#ifndef PULSAR_GUARD_SCF__BASICFOCKBUILD_HPP_
#define PULSAR_GUARD_SCF__BASICFOCKBUILD_HPP_

#include "Methods/SCF/SCFCommon.hpp"

#include <pulsar/modulebase/FockBuilder.hpp>

#include <vector>
#include <Eigen/Dense>

namespace pulsarmethods {

class BasicFockBuild : public pulsar::modulebase::FockBuilder
{
    public:
        BasicFockBuild(ID_t id) :  pulsar::modulebase::FockBuilder(id) { }
        
        virtual void initialize_(unsigned int deriv,
                                 const pulsar::datastore::Wavefunction & wfn,
                                 const pulsar::system::BasisSet & bs);

        virtual pulsar::math::IrrepSpinMatrixD calculate_(const pulsar::datastore::Wavefunction & wfn);


    private:
        std::vector<double> eri_;

        std::shared_ptr<const Eigen::MatrixXd> Hcore_;
};

}

#endif

