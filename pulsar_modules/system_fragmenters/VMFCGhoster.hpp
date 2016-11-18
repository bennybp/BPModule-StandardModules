#pragma once

#include<pulsar/modulebase/SystemFragmenter.hpp>

/** \brief Makes a ghoster consistent with a VMFC BSSE correction
 * 
 *  At the end of this you will get back your original set of fragments plus the
 *  each \f$k\f$-mer in the basis sets ranging from the \f$k+1\f$-mer basis to
 *  the \f$m\f$-mer basis, where \f$m\f$ is maximum length of a serial number in
 *  any fragment.
 */
class VMFCGhoster: public pulsar::SystemFragmenter{
public:
    ///Import the constructor from base class
        using pulsar::SystemFragmenter::SystemFragmenter;
        
        ///Actually fragments the molecule
        virtual pulsar::NMerSetType fragmentize_(const pulsar::System& mol);
};
