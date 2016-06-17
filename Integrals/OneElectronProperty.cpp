#include <pulsar/system/BasisSet.hpp>
#include <pulsar/system/AOIterator.hpp>
#include <pulsar/modulebase/OneElectronIntegral.hpp>

#include "Integrals/OneElectronProperty.hpp"


using namespace pulsar::modulemanager;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;
using namespace pulsar::modulebase;
using namespace pulsar::math;


std::vector<double>
OneElectronProperty::Calculate_(unsigned int deriv,
                                const Wavefunction & wfn,
                                const BasisSet & bs1,
                                const BasisSet & bs2)

{
    const size_t nshell1 = bs1.NShell();
    const size_t nshell2 = bs2.NShell();

    auto mod = CreateChildFromOption<OneElectronIntegral>("KEY_ONEEL_MOD");
    mod->Initialize(deriv, wfn, bs1, bs2);
    const unsigned int ncomp = mod->NComponents(); 

    // size of the workspace depends on the number of components
    size_t worksize = ncomp*bs1.MaxNFunctions()*bs2.MaxNFunctions();
    std::vector<double> work(worksize);

    std::vector<double> val(ncomp);
    std::fill(val.begin(), val.end(), 0.0);

    for(size_t n1 = 0; n1 < nshell1; n1++)
    {
        const auto & sh1 = bs1.Shell(n1);
        const size_t rowstart = bs1.ShellStart(n1);

        for(size_t n2 = 0; n2 < nshell2; n2++)
        {
            const auto & sh2  = bs2.Shell(n2);
            const size_t colstart = bs2.ShellStart(n2);

            // calculate
            size_t ncalc = mod->Calculate(n1, n2, work.data(), worksize);

            // iterate and fill in the matrix
            AOIterator<2> aoit({sh1, sh2}, false);

            // make sure the right number of integrals was returned
            if(ncalc != aoit.NFunctions())
                throw GeneralException("Bad number of integrals returned",
                                       "ncalc", ncalc, "expected", aoit.NFunctions());


            do {
                const size_t i = rowstart+aoit.ShellFunctionIdx<0>();
                const size_t j = colstart+aoit.ShellFunctionIdx<1>();

                for(auto s : wfn.opdm->GetSpins(Irrep::A))
                {
                    for(unsigned int c = 0; c < ncomp; c++)
                        val[c] += wfn.opdm->Get(Irrep::A, s)->get_value({i,j}) * work[aoit.TotalIdx() + c*ncalc];
                }
            } while(aoit.Next());
        }
    }

    return val;
}
