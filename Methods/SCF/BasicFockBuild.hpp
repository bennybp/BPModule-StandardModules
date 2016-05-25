#ifndef PULSAR_GUARD_SCF__BASICFOCKBUILD_HPP_
#define PULSAR_GUARD_SCF__BASICFOCKBUILD_HPP_

#include "SCF_Common.hpp"

#include <vector>
#include <pulsar/modulebase/FockBuilder.hpp>
#include <Eigen/Dense>

namespace pulsarmethods {

class BasicFockBuild : public pulsar::modulebase::FockBuilder
{
    public:
        BasicFockBuild(ID_t id) :  pulsar::modulebase::FockBuilder(id), initialized_(false) { }
        
        virtual pulsar::math::IrrepSpinMatrixD Build_(const pulsar::datastore::Wavefunction & wfn);

    private:
        bool initialized_;

        std::vector<double> eri_;
        Eigen::MatrixXd S12_;
        Eigen::MatrixXd Hcore_;

        void Initialize_(const pulsar::datastore::Wavefunction & wfn);
};

}

#endif

