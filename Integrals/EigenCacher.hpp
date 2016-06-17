#ifndef _GUARD_EIGENCACHER_HPP_
#define _GUARD_EIGENCACHER_HPP_

#include <pulsar/modulebase/OneElectronCacher.hpp>
#include <pulsar/system/BasisSet.hpp>

class EigenCacher : public pulsar::modulebase::OneElectronCacher
{
    public:
        using pulsar::modulebase::OneElectronCacher::OneElectronCacher;

        virtual ReturnType Calculate_(const std::string & key,
                                      unsigned int deriv,
                                      const pulsar::datastore::Wavefunction & wfn,
                                      const pulsar::system::BasisSet & bs1,
                                      const pulsar::system::BasisSet & bs2);
};


#endif
