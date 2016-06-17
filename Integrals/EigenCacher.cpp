#include "Common/EigenCommon.hpp"
#include "Integrals/EigenCacher.hpp"

using namespace pulsar::modulebase;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;


EigenCacher::ReturnType
EigenCacher::Calculate_(const std::string & key,
                        unsigned int deriv,
                        const Wavefunction & wfn,
                        const BasisSet & bs1,
                        const BasisSet & bs2)
{
    const size_t nshell1 = bs1.NShell();
    const size_t nshell2 = bs2.NShell();

    auto mod = CreateChildFromOption<OneElectronIntegral>(key);
    mod->Initialize(deriv, wfn, bs1, bs2);
    const unsigned int ncomp = mod->NComponents(); 


    std::vector<std::shared_ptr<EigenMatrixImpl>> ret(ncomp,
                                                      std::make_shared<EigenMatrixImpl>(MatrixXd(nfunc1, nfunc2)));
    
     
    return ReturnType();
}

