#ifndef PULSAR_GUARD_SCF__DIIS_HPP_
#define PULSAR_GUARD_SCF__DIIS_HPP_

#include <pulsar/modulebase/EnergyMethod.hpp>
#include <Eigen/Dense>

namespace pulsarmethods {

class DIIS : public pulsar::modulebase::EnergyMethod
{
    public:
        DIIS(ID_t id) : pulsar::modulebase::EnergyMethod(id), initialized_(false) { }
        
        virtual DerivReturnType Deriv_(size_t order, const pulsar::datastore::Wavefunction & wfn);


    private:
        bool initialized_;
        double nucrep_;
        Eigen::MatrixXd Hcore_;
        Eigen::MatrixXd S_;

        void Initialize_(const pulsar::datastore::Wavefunction & wfn);

        double CalculateEnergy_(const pulsar::math::IrrepSpinMatrixD & Dmat,
                                const pulsar::math::IrrepSpinMatrixD & Fmat);
};

}

#endif

