#ifndef _GUARD_ATOMIZER_HPP_
#define _GUARD_ATOMIZER_HPP_

#include <bpmodule/modulebase/SystemFragmenter.hpp>


/** \brief A simple fragmenter that makes each atom its own fragment
 * 
 *   RMR-Is this just a test?  Wouldn't using the iterators of the system
 *   accomplish the same thing?
 */ 
class Atomizer : public bpmodule::modulebase::SystemFragmenter
{
public:
    using SystemFragmenter::SystemFragmenter;

    virtual bpmodule::system::SystemMap Fragmentize_(const bpmodule::system::System & mol);

};

#endif
