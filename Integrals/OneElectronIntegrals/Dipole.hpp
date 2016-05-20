#ifndef _GUARD_DIPOLE_HPP_
#define _GUARD_DIPOLE_HPP_

#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include <pulsar/system/BasisSet.hpp>

class Dipole : public pulsar::modulebase::OneElectronIntegral
{
    public:
        using pulsar::modulebase::OneElectronIntegral::OneElectronIntegral;

        virtual void SetBases_(const pulsar::datastore::Wavefunction & wfn,
                               const pulsar::system::BasisSet & bs1,
                               const pulsar::system::BasisSet & bs2);

        virtual unsigned int NComponents_(void) const { return 3; }

        virtual uint64_t Calculate_(uint64_t deriv,
                                    uint64_t shell1, uint64_t shell2,
                                    double * outbuffer, size_t bufsize);

    private:
        std::vector<double> work_;

        double * transformwork_;
        double * sourcework_;
        double * xyzwork_[3];

        std::shared_ptr<pulsar::system::BasisSet> bs1_, bs2_;
};


#endif
