#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOIterator.hpp>
#include "Methods/SCF/SCF_Common.hpp"

using Eigen::MatrixXd;

using namespace pulsar::modulemanager;
using namespace pulsar::modulebase;
using namespace pulsar::system;



namespace pulsarmethods{


MatrixXd
FillOneElectronMatrix(ModulePtr<OneElectronIntegral> & mod,
                      const BasisSet & bs)
{
    size_t nao = bs.NFunctions();
    size_t nshell = bs.NShell();
    size_t maxnfunc = bs.MaxNFunctions();
    const size_t maxnfunc2 = maxnfunc * maxnfunc;

    // matrix we are returning
    MatrixXd mat(nao, nao);

    // buffer
    std::vector<double> b(maxnfunc2);

    for(size_t n1 = 0; n1 < nshell; n1++)
    {
        const auto & sh1 = bs.Shell(n1);
        const size_t rowstart = bs.ShellStart(n1);
        const size_t nfunc1 = sh1.NFunctions();

        for(size_t n2 = 0; n2 <= n1; n2++)
        {
            const auto & sh2  = bs.Shell(n2);
            const size_t colstart = bs.ShellStart(n2);
            const size_t nfunc2 = sh2.NFunctions();

            // calculate
            size_t ncalc = mod->Calculate(0, n1, n2, b.data(), maxnfunc2); 

            // iterate and fill in the matrix
            AOIterator<2> aoit({sh1, sh2}, false);

            do { 
                const size_t i = rowstart+aoit.ShellFunctionIdx<0>();
                const size_t j = colstart+aoit.ShellFunctionIdx<1>();

                mat(i,j) = mat(j, i) = b[aoit.TotalIdx()];
            } while(aoit.Next());
        }
    }

    return mat;
}


}//End namespace
