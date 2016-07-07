#pragma once

#include <pulsar/modulebase/PropertyCalculator.hpp>

namespace psr_modules {
namespace integrals {


class OneElectronProperty : public pulsar::modulebase::PropertyCalculator
{
    public:
        using pulsar::modulebase::PropertyCalculator::PropertyCalculator;

        virtual std::vector<double> calculate_(unsigned int deriv,
                                               const pulsar::datastore::Wavefunction & wfn,
                                               const pulsar::system::BasisSet & bs1,
                                               const pulsar::system::BasisSet & bs2);
};


} // close namespace integrals
} // close namespace psr_modules
