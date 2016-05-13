#include "Methods/SCF/BasicFockBuild.hpp"
#include "pulsar/modulebase/All.hpp"

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


static
BlockedEigenMatrix BuildFock(const IrrepSpinMatrixD & Dmat,
                             const MatrixXd & Hcore,
                             const std::vector<double> & eri)
{
    BlockedEigenMatrix Fmat;

    for(auto ir : Dmat.GetIrreps())
    for(auto s : Dmat.GetSpins(ir))
    {
        MappedConstMatrix D = MapConstSimpleMatrix(Dmat.Get(ir,s));

        MatrixXd F = Hcore;

        size_t nao1 = Hcore.rows();
        size_t nao2 = Hcore.cols();

        for(size_t mu = 0; mu < nao1; mu++)
        for(size_t nu = 0; nu < nao2; nu++)
        {
            for(size_t lambda = 0; lambda < nao1; lambda++)
            for(size_t sigma = 0; sigma < nao2; sigma++)
            {
                size_t mnls = INDEX4(mu, nu, lambda, sigma);
                size_t mlns = INDEX4(mu, lambda, nu, sigma);
                F(mu, nu) += 0.5*D(lambda, sigma) * (2*eri.at(mnls)-eri.at(mlns));
            }
        }

        Fmat.Take(ir, s, std::move(F));
    }

    return Fmat;
}


void BasicFockBuild::Initialize_(const Wavefunction & wfn)
{
    if(!wfn.system)
        throw GeneralException("System is not set!");

    // get the basis set
    const System & sys = *(wfn.system);
    std::string bstag = Options().Get<std::string>("BASIS_SET");

    out.Output("Obtaining basis set %? from system\n", bstag);
    const BasisSet bs = sys.GetBasisSet(bstag);

    /////////////////////////
    // Load the ERI to core
    /////////////////////////
    auto mod_ao_eri = CreateChildFromOption<TwoElectronIntegral>("KEY_AO_ERI");
    mod_ao_eri->SetBases(sys, bstag, bstag, bstag, bstag);
    eri_ = FillTwoElectronVector(mod_ao_eri, bs);


    ///////////////////////
    // Overlap
    ///////////////////////
    auto mod_ao_overlap = CreateChildFromOption<OneElectronIntegral>("KEY_AO_OVERLAP");
    mod_ao_overlap->SetBases(sys, bstag, bstag);
    MatrixXd overlap_mat = FillOneElectronMatrix(mod_ao_overlap, bs);

    // diagonalize the overlap
    SelfAdjointEigenSolver<MatrixXd> esolve(overlap_mat);
    MatrixXd s_evec = esolve.eigenvectors();
    VectorXd s_eval = esolve.eigenvalues();

    // not sure an easier way to do this
    for(size_t i = 0; i < s_eval.size(); i++)
        s_eval(i) = 1.0/sqrt(s_eval(i));

    // the S^(-1/2) matrix
    S12_ = s_evec * s_eval.asDiagonal() * s_evec.transpose();


    ////////////////////////////
    // One-electron hamiltonian
    auto mod_ao_core = CreateChildFromOption<OneElectronIntegral>("KEY_AO_COREBUILD");
    mod_ao_core->SetBases(sys, bstag, bstag);
    Hcore_ = FillOneElectronMatrix(mod_ao_core, bs);

    initialized_ = true;
}


IrrepSpinMatrixD BasicFockBuild::Build_(const Wavefunction & wfn)
{
    // has the density been set
    if(!wfn.opdm)
        throw GeneralException("Missing OPDM");

    if(!initialized_)
        Initialize_(wfn);

    // the fock matrix we are returning
    IrrepSpinMatrixD Fmat;

    for(auto ir : wfn.opdm->GetIrreps())
    for(auto s : wfn.opdm->GetSpins(ir))
    {
        MappedConstMatrix D = MapConstSimpleMatrix(wfn.opdm->Get(ir,s));

        SimpleMatrixD simpleF(D.rows(), D.cols());
        MappedMatrix F = MapSimpleMatrix(simpleF);
        F = Hcore_;

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
                F(mu, nu) += 0.5*D(lambda, sigma) * (2*eri_.at(mnls)-eri_.at(mlns));
            }
        }

        Fmat.Take(ir, s, std::move(simpleF));
    }

    return Fmat;
}


} // close namespace pulsarmethods
