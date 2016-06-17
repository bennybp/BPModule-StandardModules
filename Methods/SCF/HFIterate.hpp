#ifndef PULSAR_GUARD_SCF__HFITERATE_HPP_
#define PULSAR_GUARD_SCF__HFITERATE_HPP_

#include "Methods/SCF/SCFCommon.hpp"

#include <pulsar/modulebase/SCFIterator.hpp>
#include <Eigen/Dense>

namespace pulsarmethods {

class HFIterate : public pulsar::modulebase::SCFIterator
{
    public:
        HFIterate(ID_t id) :  pulsar::modulebase::SCFIterator(id), initialized_(false) { }
        
        virtual pulsar::datastore::Wavefunction
        Next_(const pulsar::datastore::Wavefunction & wfn, const pulsar::math::IrrepSpinMatrixD & fmat);

    private:
        bool initialized_;
        Eigen::MatrixXd S12_;

        void Initialize_(const pulsar::datastore::Wavefunction & wfn);
};

}

#endif

