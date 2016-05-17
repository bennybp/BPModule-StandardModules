#ifndef _GUARD_NUCLEARDIPOLE_HPP_
#define _GUARD_NUCLEARDIPOLE_HPP_

#include <pulsar/modulebase/SystemIntegral.hpp>

class NuclearDipole : public pulsar::modulebase::SystemIntegral
{
public:
    using pulsar::modulebase::SystemIntegral::SystemIntegral;

    virtual uint64_t Calculate_(uint64_t deriv, const pulsar::system::System & sys,
                                double * outbuffer, size_t bufsize);

};


#endif
