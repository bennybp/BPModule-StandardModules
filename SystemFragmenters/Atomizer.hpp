#ifndef _GUARD_ATOMIZER_HPP_
#define _GUARD_ATOMIZER_HPP_

#include <pulsar/modulebase/SystemFragmenter.hpp>


/** \brief A simple fragmenter that makes each atom its own fragment
 * 
 *   RMR-Is this just a test?  Wouldn't using the iterators of the system
 *   accomplish the same thing?
 */ 
class Atomizer : public pulsar::modulebase::SystemFragmenter
{
public:
    using SystemFragmenter::SystemFragmenter;

    virtual pulsar::system::SystemMap Fragmentize_(const pulsar::system::System & mol);

};

#endif
