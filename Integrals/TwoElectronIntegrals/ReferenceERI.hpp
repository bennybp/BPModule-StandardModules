#ifndef _GUARD_REFERENCERI_HPP_
#define _GUARD_REFERENCERI_HPP_

#include <pulsar/modulebase/TwoElectronIntegral.hpp>
#include <pulsar/system/BasisSet.hpp>

class ReferenceERI : public pulsar::modulebase::TwoElectronIntegral
{
public:
    ReferenceERI(ID_t id);

    virtual void SetBases_(const std::string & bs1, const std::string & bs2,
                           const std::string & bs3, const std::string & bs4);


    virtual uint64_t Calculate_(size_t deriv,
                                size_t shell1, size_t shell2,
                                size_t shell3, size_t shell4,
                                double * outbuffer, size_t bufsize);

    virtual ~ReferenceERI();


private:
    std::shared_ptr<pulsar::system::BasisSet> bs1_, bs2_, bs3_, bs4_;

};


#endif
