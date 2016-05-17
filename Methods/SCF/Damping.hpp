#ifndef PULSAR_GUARD_SCF__BPTEST_HPP_
#define PULSAR_GUARD_SCF__BPTEST_HPP_

#include <pulsar/modulebase/EnergyMethod.hpp>
#include <eigen3/Eigen/Dense>

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

        void Initialize_(const pulsar::system::System & sys, const std::string & bstag);

        double CalculateEnergy_(const pulsar::math::IrrepSpinMatrixD & Dmat,
                                const pulsar::math::IrrepSpinMatrixD & Fmat);
};

}

#endif

