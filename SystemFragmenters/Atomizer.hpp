#ifndef _GUARD_ATOMIZER_HPP_
#define _GUARD_ATOMIZER_HPP_

#include <pulsar/modulebase/SystemFragmenter.hpp>


/** \brief A simple fragmenter that makes each atom its own fragment
 * 
 */ 
class Atomizer : public pulsar::modulebase::SystemFragmenter
{
public:
    using SystemFragmenter::SystemFragmenter;

    virtual NMerSetType fragmentize_(const pulsar::system::System & mol);

};

#endif
