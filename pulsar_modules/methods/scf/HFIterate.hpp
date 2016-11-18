#ifndef PULSAR_GUARD_SCF__HFITERATE_HPP_
#define PULSAR_GUARD_SCF__HFITERATE_HPP_

#include "Methods/SCF/SCFCommon.hpp"

#include <pulsar/modulebase/SCFIterator.hpp>
#include <Eigen/Dense>

namespace pulsarmethods {

class HFIterate : public pulsar::SCFIterator
{
    public:
        HFIterate(ID_t id) :  pulsar::SCFIterator(id), initialized_(false) { }
        
        virtual pulsar::Wavefunction
        next_(const pulsar::Wavefunction & wfn, const pulsar::IrrepSpinMatrixD & fmat);

    private:
        bool initialized_;
        Eigen::MatrixXd S12_;

        void initialize_(const pulsar::Wavefunction & wfn);
};

}

#endif

