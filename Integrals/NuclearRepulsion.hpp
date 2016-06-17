#ifndef _GUARD_NUCLEARREPULSION_HPP_
#define _GUARD_NUCLEARREPULSION_HPP_

#include <pulsar/modulebase/SystemIntegral.hpp>

class NuclearRepulsion : public pulsar::modulebase::SystemIntegral
{
    public:
        using pulsar::modulebase::SystemIntegral::SystemIntegral;

        virtual void Initialize_(unsigned int deriv, const pulsar::system::System & sys);

        virtual uint64_t Calculate_(double * outbuffer, size_t bufsize);

    private:
        const pulsar::system::System * sys_;
};


#endif
