#include "Integrals/OneElectronIntegralSum.hpp"

#include <pulsar/util/StringUtil.hpp>

using namespace pulsar::exception;
using namespace pulsar::datastore;
using namespace pulsar::system;


namespace psr_modules {
namespace integrals {


uint64_t OneElectronIntegralSum::calculate_(uint64_t shell1, uint64_t shell2,
                                            double * outbuffer, size_t bufsize)
{
    // Peel the first one so we can get a guess for the number of elements.
    // We can place this one directly into the output buffer
    auto it = modules_.begin();

    uint64_t n_initial = it->second->calculate(shell1, shell2,
                                               outbuffer, bufsize);

    // now allocate the buffer and loop over the rest
    std::vector<double> tmpbuf(n_initial, 0.0);

    // loop over the rest
    ++it;
    while(it != modules_.end())
    {
        uint64_t n = it->second->calculate(shell1, shell2,
                                           tmpbuf.data(),
                                           n_initial);

        if(n != n_initial)
            throw GeneralException("Error - inconsistent number of values returned by OneElectronIntegrals",
                                   "n", n, "nexpected", n_initial,
                                   "modulekey", it->second->key(),
                                   "modulename", it->second->name());

        for(size_t i = 0; i < n; i++)
            outbuffer[i] += tmpbuf[i];

        ++it;
    }

    return n_initial;
}



void OneElectronIntegralSum::initialize_(unsigned int deriv,
                                         const Wavefunction & wfn,
                                         const BasisSet & bs1,
                                         const BasisSet & bs2)
{
    using pulsar::util::line;

    const auto mods = options().get<std::vector<std::string>>("KEY_AO_TERMS");
    if(mods.size() == 0)
        throw GeneralException("No modules given to OneElectronIntegralSum");


    for(const auto & a : mods)
    {
        auto mod_add = create_child<OneElectronIntegral>(a);
        mod_add->initialize(deriv, wfn, bs1, bs2);
        modules_.emplace(a, std::move(mod_add));
    }


    // Print out my info
    out.output("%? initialized with %? modules\n",
               this->name(), modules_.size());

    out.output("    %-20?  %4?  %-20?\n", "Module Key", "ID", "Module Name");
    out.output(line('-'));
    for(const auto & mod : modules_)
        out.output("    %-20?  %4?  %-20?\n", mod.second->key(), mod.second->id(), mod.second->name());
    out.output(line('-'));
    out.output("\n");
}

} // close namespace integrals
} // close namespace psr_modules
