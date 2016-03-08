#ifndef _GUARD_ATOMIZER_HPP_
#define _GUARD_ATOMIZER_HPP_

#include <bpmodule/modulebase/SystemFragmenter.hpp>

class Atomizer : public bpmodule::modulebase::SystemFragmenter
{
public:
    using SystemFragmenter::SystemFragmenter;

    virtual bpmodule::system::SystemMap Fragmentize_(const bpmodule::system::System & mol);

};

#endif
