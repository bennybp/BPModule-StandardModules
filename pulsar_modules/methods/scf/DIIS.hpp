#ifndef PULSAR_GUARD_SCF__DIIS_HPP_
#define PULSAR_GUARD_SCF__DIIS_HPP_

#include <pulsar/modulebase/EnergyMethod.hpp>
#include <Eigen/Dense>

namespace pulsarmethods {

class DIIS : public pulsar::EnergyMethod
{
    public:
        DIIS(ID_t id) : pulsar::EnergyMethod(id), initialized_(false) { }
        
        virtual pulsar::DerivReturnType deriv_(size_t order, const pulsar::Wavefunction & wfn);


    private:
        bool initialized_;
        double nucrep_;
        std::shared_ptr<const Eigen::MatrixXd> Hcore_;
        std::shared_ptr<const Eigen::MatrixXd> S_;

        void initialize_(const pulsar::Wavefunction & wfn);

        double CalculateEnergy_(const pulsar::IrrepSpinMatrixD & Dmat,
                                const pulsar::IrrepSpinMatrixD & Fmat);
};

}

#endif

