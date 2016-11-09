#pragma once

#include <pulsar/modulebase/SystemFragmenter.hpp>


/** \brief A simple fragmenter that makes each atom its own fragment
 * 
 */ 
class Atomizer : public pulsar::SystemFragmenter
{
public:
    using SystemFragmenter::SystemFragmenter;

    virtual pulsar::NMerSetType fragmentize_(const pulsar::System & mol);

};
