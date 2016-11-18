#ifndef PULSAR_GUARD_SCF__BASICFOCKBUILD_HPP_
#define PULSAR_GUARD_SCF__BASICFOCKBUILD_HPP_

#include "Methods/SCF/SCFCommon.hpp"

#include <pulsar/modulebase/FockBuilder.hpp>

#include <vector>
#include <Eigen/Dense>

namespace pulsarmethods {

class BasicFockBuild : public pulsar::FockBuilder
{
    public:
        BasicFockBuild(ID_t id) :  pulsar::FockBuilder(id) { }
        
        virtual void initialize_(unsigned int deriv,
                                 const pulsar::Wavefunction & wfn,
                                 const pulsar::BasisSet & bs);

        virtual pulsar::IrrepSpinMatrixD calculate_(const pulsar::Wavefunction & wfn);


    private:
        std::vector<double> eri_;

        std::shared_ptr<const Eigen::MatrixXd> Hcore_;
};

}

#endif

