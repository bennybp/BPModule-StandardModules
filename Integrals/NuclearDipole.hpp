#ifndef _GUARD_NUCLEARDIPOLE_HPP_
#define _GUARD_NUCLEARDIPOLE_HPP_

#include <pulsar/modulebase/SystemIntegral.hpp>

class NuclearDipole : public pulsar::modulebase::SystemIntegral
{
    public:
        using pulsar::modulebase::SystemIntegral::SystemIntegral;

        virtual void initialize_(unsigned int deriv, const pulsar::system::System & sys);

        virtual uint64_t calculate_(double * outbuffer, size_t bufsize);

    private:
        const pulsar::system::System * sys_;
};


#endif
