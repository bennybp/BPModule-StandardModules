#include "Methods/SCF/BasicFockBuild.hpp"

#include <pulsar/modulebase/All.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;

using namespace pulsar::datastore;
using namespace pulsar::system;
using namespace pulsar::math;
using namespace pulsar::exception;
using namespace pulsar::output;
using namespace pulsar::modulebase;


namespace pulsarmethods {


void BasicFockBuild::Initialize_(unsigned int deriv, const Wavefunction & wfn,
                                 const BasisSet & bs)
{
    if(!wfn.system)
        throw GeneralException("System is not set!");

    /////////////////////////
    // Load the ERI to core
    /////////////////////////
    auto mod_ao_eri = CreateChildFromOption<TwoElectronIntegral>("KEY_AO_ERI");
    mod_ao_eri->Initialize(0, wfn, bs, bs, bs, bs);
    eri_ = FillTwoElectronVector(mod_ao_eri, bs);


    ///////////////////////
    // Overlap
    ///////////////////////
    auto mod_ao_overlap = CreateChildFromOption<OneElectronIntegral>("KEY_AO_OVERLAP");
    mod_ao_overlap->Initialize(0, wfn, bs, bs);
    MatrixXd overlap_mat = FillOneElectronMatrix(mod_ao_overlap, bs);

    // diagonalize the overlap
    SelfAdjointEigenSolver<MatrixXd> esolve(overlap_mat);
    MatrixXd s_evec = esolve.eigenvectors();
    VectorXd s_eval = esolve.eigenvalues();

    // not sure an easier way to do this
    for(int i = 0; i < s_eval.size(); i++)
        s_eval(i) = 1.0/sqrt(s_eval(i));

    // the S^(-1/2) matrix
    S12_ = s_evec * s_eval.asDiagonal() * s_evec.transpose();


    ////////////////////////////
    // One-electron hamiltonian
    ///////////////////////
    auto mod_ao_core = CreateChildFromOption<OneElectronIntegral>("KEY_AO_COREBUILD");
    mod_ao_core->Initialize(0, wfn, bs, bs);
    Hcore_ = FillOneElectronMatrix(mod_ao_core, bs);
}


IrrepSpinMatrixD BasicFockBuild::Calculate_(const Wavefunction & wfn)
{
    if(!wfn.opdm)
        throw GeneralException("Missing OPDM");

    const size_t nao1 = Hcore_.rows();
    const size_t nao2 = Hcore_.cols();


    // the fock matrix we are returning
    IrrepSpinMatrixD Fmat;

    for(auto ir : wfn.opdm->GetIrreps())
    {
        const auto & spins = wfn.opdm->GetSpins(ir);
        if(spins == std::set<int>{0})
        {
            // Restricted
            std::shared_ptr<const MatrixXd> Dptr = convert_to_eigen(wfn.opdm->Get(ir,0));
            const auto & D = *Dptr;

            MatrixXd F(Hcore_);

            for(size_t mu = 0; mu < nao1; mu++)
            for(size_t nu = 0; nu < nao2; nu++)
            {
                for(size_t lambda = 0; lambda < nao1; lambda++)
                for(size_t sigma = 0; sigma < nao2; sigma++)
                {
                    size_t mnls = INDEX4(mu, nu, lambda, sigma);
                    size_t mlns = INDEX4(mu, lambda, nu, sigma);
                    F(mu, nu) += 0.5*D(lambda, sigma) * (2*eri_.at(mnls)-eri_.at(mlns));
                }
            }

            auto Fimpl = std::make_shared<EigenMatrixImpl>(std::move(F));
            Fmat.Take(ir, 0, std::move(Fimpl));
        }
        else
        {
            /*
            // unrestricted
            MappedConstMatrix Dalpha = MapConstSimpleMatrix(wfn.opdm->Get(ir, 1));
            MappedConstMatrix Dbeta = MapConstSimpleMatrix(wfn.opdm->Get(ir, -1));

            SimpleMatrixD simpleFalpha(Dalpha.rows(), Dalpha.cols());
            SimpleMatrixD simpleFbeta(Dbeta.rows(), Dbeta.cols());
            MappedMatrix Falpha = MapSimpleMatrix(simpleFalpha);
            MappedMatrix Fbeta = MapSimpleMatrix(simpleFbeta);
            Falpha = Hcore_;
            Fbeta = Hcore_;

            size_t nao1 = Hcore_.rows();
            size_t nao2 = Hcore_.cols();

            for(size_t mu = 0; mu < nao1; mu++)
            for(size_t nu = 0; nu < nao2; nu++)
            {
                for(size_t lambda = 0; lambda < nao1; lambda++)
                for(size_t sigma = 0; sigma < nao2; sigma++)
                {
                    size_t mnls = INDEX4(mu, nu, lambda, sigma);
                    size_t mlns = INDEX4(mu, lambda, nu, sigma);
                    Falpha(mu, nu) += (Dalpha(lambda, sigma) * (eri_.at(mnls)-eri_.at(mlns)) + Dbeta(lambda, sigma)*eri_.at(mnls));
                    Fbeta(mu, nu) += (Dbeta(lambda, sigma) * (eri_.at(mnls)-eri_.at(mlns)) + Dalpha(lambda, sigma)*eri_.at(mnls));
                }
            }

            Fmat.Take(ir, 1, std::move(simpleFalpha));
            Fmat.Take(ir, -1, std::move(simpleFbeta));
            */
        }
    }

    return Fmat;
}


} // close namespace pulsarmethods
