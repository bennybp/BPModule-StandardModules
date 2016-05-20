#ifndef _GUARD_ONEELECTRONPROPERTY_HPP_
#define _GUARD_ONEELECTRONPROPERTY_HPP_

#include <pulsar/modulebase/PropertyCalculator.hpp>

class OneElectronProperty : public pulsar::modulebase::PropertyCalculator
{
    public:
        using pulsar::modulebase::PropertyCalculator::PropertyCalculator;

        virtual std::vector<double> Calculate_(const pulsar::datastore::Wavefunction & wfn,
                                               const pulsar::system::BasisSet & bs1,
                                               const pulsar::system::BasisSet & bs2);
};


#endif
