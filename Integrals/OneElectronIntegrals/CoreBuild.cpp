#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/NShellFunction.hpp>
#include <pulsar/constants.h>

#include "../Common.hpp"
#include "CoreBuild.hpp"


// Get a value of S_IJ
#define S_IJ(i,j) (s_ij[((i)*(nam2) + j)])

using namespace pulsar::modulemanager;
using namespace pulsar::modulebase;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;


static void CallAndBuild(ModulePtr<OneElectronIntegral> & mod,
                         double * tmpbuf, size_t nexpected,
                         uint64_t shell1, uint64_t shell2,
                         double * outbuffer, size_t bufsize)
{
    uint64_t n = mod->Calculate(shell1, shell2,
                                tmpbuf, bufsize);


    if(n != nexpected)
        throw GeneralException("Error - inconsistent number of values returned by OneElectronIntegrals",
                               "n", n, "nexpected", nexpected,
                               "modulekey", mod->Key(), "modulename", mod->Name());

    for(size_t i = 0; i < n; i++)
        outbuffer[i] += tmpbuf[i];
}


uint64_t CoreBuild::Calculate_(uint64_t shell1, uint64_t shell2,
                               double * outbuffer, size_t bufsize)
{
    if(modules_.size() == 0)
        throw GeneralException("I don't have any modules?");

    auto it = modules_.begin();

    // put the first directly in the output buffer
    uint64_t n_initial = it->second->Calculate(shell1, shell2,
                                               outbuffer, bufsize);

    // now allocate the buffer and loop over the rest
    std::vector<double> tmpbuf(n_initial, 0.0);

    // loop over the rest
    ++it;
    while(it != modules_.end())
    {
        CallAndBuild(it->second,
                     tmpbuf.data(), n_initial,
                     shell1, shell2, outbuffer, bufsize);
        ++it;
    }

    return n_initial;
}



void CoreBuild::Initialize_(unsigned int deriv,
                            const Wavefunction & wfn,
                            const BasisSet & bs1,
                            const BasisSet & bs2)
{
    ///////////////////////////////// 
    // load all the required modules
    ///////////////////////////////// 

    /////////////////////// 
    // Kinetic Energy
    auto mod_ao_kinetic = CreateChildFromOption<OneElectronIntegral>("KEY_AO_KINETIC");
    mod_ao_kinetic->Initialize(deriv, wfn, bs1, bs2);
    modules_.emplace("Kinetic Energy", std::move(mod_ao_kinetic));

    /////////////////////// 
    // Nuclear Attraction
    auto mod_ao_nucatt = CreateChildFromOption<OneElectronIntegral>("KEY_AO_NUCATT");
    mod_ao_nucatt->Initialize(deriv, wfn, bs1, bs2);
    modules_.emplace("Electron-Nuclear Attraction", std::move(mod_ao_nucatt));

    // do the additional terms
    if(Options().Has("KEY_AO_ADDITIONAL"))
    {
        const auto add = Options().Get<std::vector<std::string>>("KEY_AO_ADDITIONAL");
        for(const auto & a : add)
        {
            auto mod_add = CreateChild<OneElectronIntegral>(a);
            mod_add->Initialize(deriv, wfn, bs1, bs2);
            modules_.emplace(a, std::move(mod_add));
        }
    }
}
