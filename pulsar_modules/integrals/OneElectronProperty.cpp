#include <pulsar/system/BasisSet.hpp>
#include <pulsar/system/AOIterator.hpp>
#include <pulsar/modulebase/OneElectronIntegral.hpp>

#include "Integrals/OneElectronProperty.hpp"


using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;
using namespace pulsar::modulebase;
using namespace pulsar::math;


namespace psr_modules {
namespace integrals {


std::vector<double>
OneElectronProperty::calculate_(unsigned int deriv,
                                const Wavefunction & wfn,
                                const BasisSet & bs1,
                                const BasisSet & bs2)

{
    const size_t nshell1 = bs1.n_shell();
    const size_t nshell2 = bs2.n_shell();

    auto mod = create_child_from_option<OneElectronIntegral>("KEY_ONEEL_MOD");
    mod->initialize(deriv, wfn, bs1, bs2);
    const unsigned int ncomp = mod->n_components(); 

    // size of the workspace depends on the number of components
    size_t worksize = ncomp*bs1.max_n_functions()*bs2.max_n_functions();
    std::vector<double> work(worksize);

    std::vector<double> val(ncomp);
    std::fill(val.begin(), val.end(), 0.0);

    for(size_t n1 = 0; n1 < nshell1; n1++)
    {
        const auto & sh1 = bs1.shell(n1);
        const size_t rowstart = bs1.shell_start(n1);

        for(size_t n2 = 0; n2 < nshell2; n2++)
        {
            const auto & sh2  = bs2.shell(n2);
            const size_t colstart = bs2.shell_start(n2);

            // calculate
            size_t ncalc = mod->calculate(n1, n2, work.data(), worksize);

            // iterate and fill in the matrix
            AOIterator<2> aoit({sh1, sh2}, false);

            // make sure the right number of integrals was returned
            if(ncalc != aoit.n_functions())
                throw PulsarException("Bad number of integrals returned",
                                       "ncalc", ncalc, "expected", aoit.n_functions());


            do {
                const size_t i = rowstart+aoit.shell_function_idx<0>();
                const size_t j = colstart+aoit.shell_function_idx<1>();

                for(auto s : wfn.opdm->get_spins(Irrep::A))
                {
                    for(unsigned int c = 0; c < ncomp; c++)
                        val[c] += wfn.opdm->get(Irrep::A, s)->get_value({i,j}) * work[aoit.total_idx() + c*ncalc];
                }
            } while(aoit.next());
        }
    }

    return val;
}


} // close namespace integrals
} // close namespace psr_modules
