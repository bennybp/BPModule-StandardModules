#pragma once

#include <pulsar/modulebase/OneElectronMatrix.hpp>

namespace psr_modules {
namespace integrals {


class OneElectron_Eigen : public pulsar::modulebase::OneElectronMatrix
{
    public:
        using pulsar::modulebase::OneElectronMatrix::OneElectronMatrix;

        virtual ReturnType calculate_(const std::string & key,
                                      unsigned int deriv,
                                      const pulsar::datastore::Wavefunction & wfn,
                                      const pulsar::system::BasisSet & bs1,
                                      const pulsar::system::BasisSet & bs2);
};


} // close namespace integrals
} // close namespace psr_modules
