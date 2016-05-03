#ifndef _GUARD_ONEELECTRONPOTENTIAL_HPP_
#define _GUARD_ONEELECTRONPOTENTIAL_HPP_

#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/math/Grid.hpp>



class OneElectronPotential : public pulsar::modulebase::OneElectronIntegral
{
    public:
        OneElectronPotential(ID_t id);

        virtual void SetBases_(const std::string & bs1, const std::string & bs2);

        virtual uint64_t Calculate_(uint64_t deriv, uint64_t shell1, uint64_t shell2,
                                    double * outbuffer, size_t bufsize);

        virtual ~OneElectronPotential();

    private:
        std::vector<double> work_;

        // amwork_[i][j] = work for am pair i,j
        std::vector<std::vector<double *>> amwork_;

        double * transformwork_;
        double * sourcework_;

        std::shared_ptr<pulsar::system::BasisSet> bs1_, bs2_;

        uint64_t CalculateWithGrid_(uint64_t deriv,
                                    uint64_t shell1, uint64_t shell2,
                                    const pulsar::math::Grid & grid,
                                    double * outbuffer, size_t bufsize);
};


#endif
