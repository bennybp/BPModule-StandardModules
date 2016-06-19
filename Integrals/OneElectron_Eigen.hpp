#ifndef _GUARD_ONEELECTRON_EIGEN_HPP_
#define _GUARD_ONEELECTRON_EIGEN_HPP_

#include <pulsar/modulebase/OneElectronMatrix.hpp>

class OneElectron_Eigen : public pulsar::modulebase::OneElectronMatrix
{
    public:
        using pulsar::modulebase::OneElectronMatrix::OneElectronMatrix;

        virtual ReturnType Calculate_(const std::string & key,
                                      unsigned int deriv,
                                      const pulsar::datastore::Wavefunction & wfn,
                                      const pulsar::system::BasisSet & bs1,
                                      const pulsar::system::BasisSet & bs2);
};


#endif
