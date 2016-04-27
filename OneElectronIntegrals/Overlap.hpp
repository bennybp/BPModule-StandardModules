#ifndef _GUARD_OVERLAP_HPP_
#define _GUARD_OVERLAP_HPP_

#include <pulsar/modulebase/OneElectronIntegral.hpp>

class Overlap : public pulsar::modulebase::OneElectronIntegral
{
    public:
        Overlap(ID_t id);

        virtual void SetBases_(const std::string & bs1, const std::string & bs2);

        virtual uint64_t Calculate_(uint64_t deriv, uint64_t shell1, uint64_t shell2,
                                    double * outbuffer, size_t bufsize);

        virtual ~Overlap();

    private:
        //! \todo combine all these into a single memory allocation
        std::unique_ptr<double[]> work_;
        std::unique_ptr<double[]> transformwork_;
        std::unique_ptr<double[]> srcwork_;
        double * work3_[3];

        std::shared_ptr<pulsar::system::BasisSet> bs1_, bs2_;

        std::shared_ptr<pulsar::system::BasisSet> Normalize_(const pulsar::system::BasisSet & bs);
};


#endif
