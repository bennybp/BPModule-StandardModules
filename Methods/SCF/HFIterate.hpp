#ifndef HFITERATE_HPP
#define HFITERATE_HPP

#include "SCF_Common.hpp"

#include <vector>
#include <pulsar/modulebase/EnergyMethod.hpp>
#include <eigen3/Eigen/Dense>

namespace pulsarmethods {

class HFIterate : public pulsar::modulebase::EnergyMethod
{
    public:
        HFIterate(ID_t id) :  pulsar::modulebase::EnergyMethod(id), initialized_(false) { }
        
        virtual DerivReturnType Deriv_(size_t order, const pulsar::datastore::Wavefunction & wfn);

    private:
        bool initialized_;

        std::vector<double> eri_;
        double nucrep_;
        Eigen::MatrixXd S12_;
        Eigen::MatrixXd Hcore_;


        void Initialize_(const pulsar::system::System & sys, const std::string & bstag);

        double CalculateEnergy(const pulsar::math::IrrepSpinMatrixD & Dmat,
                               const BlockedEigenMatrix & Fmat);
};

}

#endif

