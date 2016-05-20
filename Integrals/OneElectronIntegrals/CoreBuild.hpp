#ifndef _GUARD_COREBUILD_HPP_
#define _GUARD_COREBUILD_HPP_

#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include <pulsar/system/BasisSet.hpp>

class CoreBuild : public pulsar::modulebase::OneElectronIntegral
{
    public:
        using pulsar::modulebase::OneElectronIntegral::OneElectronIntegral;

        virtual void SetBases_(const pulsar::datastore::Wavefunction & wfn,
                               const pulsar::system::BasisSet & bs1,
                               const pulsar::system::BasisSet & bs2);

        virtual uint64_t Calculate_(uint64_t deriv, uint64_t shell1, uint64_t shell2,
                                    double * outbuffer, size_t bufsize);

    private:
        typedef pulsar::modulemanager::ModulePtr<pulsar::modulebase::OneElectronIntegral> OneInt;

        std::map<std::string, OneInt> modules_;
};


#endif
