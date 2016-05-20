#ifndef _GUARD_ONEELECTRONPOTENTIAL_HPP_
#define _GUARD_ONEELECTRONPOTENTIAL_HPP_

#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/math/Grid.hpp>



class OneElectronPotential : public pulsar::modulebase::OneElectronIntegral
{
    public:
        using pulsar::modulebase::OneElectronIntegral::OneElectronIntegral;

        virtual void SetBases_(const pulsar::datastore::Wavefunction & wfn,
                               const pulsar::system::BasisSet & bs1,
                               const pulsar::system::BasisSet & bs2);

        virtual uint64_t Calculate_(uint64_t deriv, uint64_t shell1, uint64_t shell2,
                                    double * outbuffer, size_t bufsize);

    private:
        std::vector<double> work_;

        // amwork_[i][j] = work for am pair i,j
        std::vector<std::vector<double *>> amwork_;
        std::shared_ptr<const pulsar::system::System> sys_;

        double * transformwork_;
        double * sourcework_;

        std::shared_ptr<pulsar::system::BasisSet> bs1_, bs2_;

        uint64_t CalculateWithGrid_(uint64_t deriv,
                                    uint64_t shell1, uint64_t shell2,
                                    const pulsar::math::Grid & grid,
                                    double * outbuffer, size_t bufsize);
};


#endif
