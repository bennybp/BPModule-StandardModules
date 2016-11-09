#pragma once

#include<pulsar/modulebase/SystemFragmenter.hpp>

/** \brief Makes a ghoster consistent with a typical CP correction
 * 
 *  At the end of this each system will be in the supersystem basis set.
 */
class CPGhoster: public pulsar::SystemFragmenter{
public:
    ///Import the constructor from base class
        using pulsar::SystemFragmenter::SystemFragmenter;
        
        ///Actually fragments the molecule
        virtual pulsar::NMerSetType fragmentize_(const pulsar::System& mol);
};
