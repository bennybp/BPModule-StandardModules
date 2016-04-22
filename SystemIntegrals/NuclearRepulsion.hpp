#ifndef _GUARD_NUCLEARREPULSION_HPP_
#define _GUARD_NUCLEARREPULSION_HPP_

#include <pulsar/modulebase/SystemIntegral.hpp>

class NuclearRepulsion : public pulsar::modulebase::SystemIntegral
{
public:
    NuclearRepulsion(ID_t id);

    virtual uint64_t Calculate_(uint64_t deriv, double * outbuffer, size_t bufsize);

    virtual ~NuclearRepulsion();

};


#endif
