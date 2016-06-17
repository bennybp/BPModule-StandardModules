#ifndef PULSAR_GUARD_SCF__BPTEST_HPP_
#define PULSAR_GUARD_SCF__BPTEST_HPP_

#include <pulsar/modulebase/EnergyMethod.hpp>
#include <Eigen/Dense>

namespace pulsarmethods {

class Damping : public pulsar::modulebase::EnergyMethod
{
    public:
        Damping(ID_t id) : pulsar::modulebase::EnergyMethod(id), initialized_(false) { }
        
        virtual DerivReturnType Deriv_(size_t order, const pulsar::datastore::Wavefunction & wfn);


    private:
        bool initialized_;
        double nucrep_;
        Eigen::MatrixXd Hcore_;

        void Initialize_(const pulsar::datastore::Wavefunction & wfn);
};

}

#endif

