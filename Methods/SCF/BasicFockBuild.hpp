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
        BasicFockBuild(ID_t id) :  pulsar::modulebase::FockBuilder(id) { }
        
        virtual void Initialize_(unsigned int deriv,
                                 const pulsar::datastore::Wavefunction & wfn,
                                 const pulsar::system::BasisSet & bs);

        virtual pulsar::math::IrrepSpinMatrixD Calculate_(const pulsar::datastore::Wavefunction & wfn);


    private:
        std::vector<double> eri_;
        Eigen::MatrixXd S12_;
        Eigen::MatrixXd Hcore_;
};

}

#endif

