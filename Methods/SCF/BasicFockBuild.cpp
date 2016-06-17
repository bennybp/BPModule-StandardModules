#include "Methods/SCF/BasicFockBuild.hpp"

#include <pulsar/modulebase/All.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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


    /////////////////////////////////////
    // The one-electron integral cacher
    /////////////////////////////////////
    auto mod_ao_cache = CreateChildFromOption<OneElectronCacher>("KEY_AO_CACHER");

    ///////////////////////
    // Overlap
    ///////////////////////
    const std::string ao_overlap_key = Options().Get<std::string>("KEY_AO_OVERLAP");
    auto overlapimpl = mod_ao_cache->Calculate(ao_overlap_key, 0, wfn, bs, bs);
    std::shared_ptr<const MatrixXd> overlap_mat = convert_to_eigen(overlapimpl.at(0));  // .at(0) = first (and only) component
    S12_ = FormS12(*overlap_mat);

    ////////////////////////////
    // One-electron hamiltonian
    ///////////////////////
    const std::string ao_build_key = Options().Get<std::string>("KEY_AO_COREBUILD");
    auto Hcoreimpl = mod_ao_cache->Calculate(ao_build_key, 0, wfn, bs, bs);
    Hcore_ = convert_to_eigen(Hcoreimpl.at(0));  // .at(0) = first (and only) component
}


IrrepSpinMatrixD BasicFockBuild::Calculate_(const Wavefunction & wfn)
{
    if(!wfn.opdm)
        throw GeneralException("Missing OPDM");

    const size_t nao1 = Hcore_->rows();
    const size_t nao2 = Hcore_->cols();


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

            MatrixXd F(*Hcore_);

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
