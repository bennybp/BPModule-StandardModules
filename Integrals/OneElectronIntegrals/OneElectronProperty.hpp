#ifndef _GUARD_ONEELECTRONPROPERTY_HPP_
#define _GUARD_ONEELECTRONPROPERTY_HPP_

#include <pulsar/modulebase/PropertyCalculator.hpp>

class OneElectronProperty : public pulsar::modulebase::PropertyCalculator
{
    public:
        using pulsar::modulebase::PropertyCalculator::PropertyCalculator;

        virtual std::vector<double> Calculate_(const pulsar::datastore::Wavefunction & wfn,
                                               const std::string & bs1,
                                               const std::string & bs2);
};


#endif
