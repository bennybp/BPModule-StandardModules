#include "Methods/SCF/BasicFockBuild.hpp"

#include <pulsar/modulebase/All.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace pulsar;

namespace pulsarmethods {


void BasicFockBuild::initialize_(unsigned int deriv, const Wavefunction & wfn,
                                 const BasisSet & bs)
{
    if(!wfn.system)
        throw PulsarException("System is not set!");

    /////////////////////////
    // Load the ERI to core
    /////////////////////////
    auto mod_ao_eri = create_child_from_option<TwoElectronIntegral>("KEY_AO_ERI");
    mod_ao_eri->initialize(0, wfn, bs, bs, bs, bs);
    eri_ = FillTwoElectronVector(mod_ao_eri, bs);


    /////////////////////////////////////
    // The one-electron integral cacher
    /////////////////////////////////////
    auto mod_ao_cache = create_child_from_option<OneElectronMatrix>("KEY_ONEEL_MAT");

    ////////////////////////////
    // One-electron hamiltonian
    ///////////////////////
    const std::string ao_build_key = options().get<std::string>("KEY_AO_COREBUILD");
    auto Hcoreimpl = mod_ao_cache->calculate(ao_build_key, 0, wfn, bs, bs);
    Hcore_ = convert_to_eigen(Hcoreimpl.at(0));  // .at(0) = first (and only) component
}


IrrepSpinMatrixD BasicFockBuild::calculate_(const Wavefunction & wfn)
{
    if(!wfn.opdm)
        throw PulsarException("Missing OPDM");

    const size_t nao1 = Hcore_->rows();
    const size_t nao2 = Hcore_->cols();


    // the fock matrix we are returning
    IrrepSpinMatrixD Fmat;

    for(auto ir : wfn.opdm->get_irreps())
    {
        const auto & spins = wfn.opdm->get_spins(ir);
        if(spins == std::set<int>{0})
        {
            // Restricted
            std::shared_ptr<const MatrixXd> Dptr = convert_to_eigen(wfn.opdm->get(ir,0));
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
            Fmat.set(ir, 0, std::move(Fimpl));
        }
        else
        {
            /*
            // unrestricted
            MappedConstMatrix Dalpha = MapConstSimpleMatrix(wfn.opdm->get(ir, 1));
            MappedConstMatrix Dbeta = MapConstSimpleMatrix(wfn.opdm->get(ir, -1));

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

            Fmat.set(ir, 1, std::move(simpleFalpha));
            Fmat.set(ir, -1, std::move(simpleFbeta));
            */
        }
    }

    return Fmat;
}


} // close namespace pulsarmethods
