#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOOrdering.hpp>
#include <pulsar/system/NShellFunction.hpp>
#include <pulsar/constants.h>

#include "Common/BasisSetCommon.hpp"
#include "Integrals/CoreBuild.hpp"


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
    uint64_t n = mod->calculate(shell1, shell2,
                                tmpbuf, bufsize);


    if(n != nexpected)
        throw GeneralException("Error - inconsistent number of values returned by OneElectronIntegrals",
                               "n", n, "nexpected", nexpected,
                               "modulekey", mod->key(), "modulename", mod->name());

    for(size_t i = 0; i < n; i++)
        outbuffer[i] += tmpbuf[i];
}


uint64_t CoreBuild::calculate_(uint64_t shell1, uint64_t shell2,
                               double * outbuffer, size_t bufsize)
{
    if(modules_.size() == 0)
        throw GeneralException("I don't have any modules?");

    auto it = modules_.begin();

    // put the first directly in the output buffer
    uint64_t n_initial = it->second->calculate(shell1, shell2,
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



void CoreBuild::initialize_(unsigned int deriv,
                            const Wavefunction & wfn,
                            const BasisSet & bs1,
                            const BasisSet & bs2)
{
    const auto add = options().get<std::vector<std::string>>("KEY_AO_CORE_TERMS");
    for(const auto & a : add)
    {
        auto mod_add = create_child<OneElectronIntegral>(a);
        mod_add->initialize(deriv, wfn, bs1, bs2);
        modules_.emplace(a, std::move(mod_add));
    }
}
