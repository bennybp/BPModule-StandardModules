#ifndef _GUARD_NUCLEARREPULSION_HPP_
#define _GUARD_NUCLEARREPULSION_HPP_

#include <pulsar/modulebase/SystemIntegral.hpp>

class NuclearRepulsion : public pulsar::modulebase::SystemIntegral
{
public:
    using pulsar::modulebase::SystemIntegral::SystemIntegral;

    virtual uint64_t Calculate_(uint64_t deriv, const pulsar::system::System & sys,
                                double * outbuffer, size_t bufsize);

};


#endif
