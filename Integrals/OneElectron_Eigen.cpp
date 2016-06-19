#include <pulsar/system/AOIterator.hpp>
#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include "Common/EigenCommon.hpp"
#include "Integrals/OneElectron_Eigen.hpp"

using Eigen::MatrixXd;

using namespace pulsar::modulebase;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::datastore;
using namespace bphash;


OneElectron_Eigen::ReturnType
OneElectron_Eigen::Calculate_(const std::string & key,
                              unsigned int deriv,
                              const Wavefunction & wfn,
                              const BasisSet & bs1,
                              const BasisSet & bs2)
{
    const bool usecache = Options().Get<bool>("CACHE_RESULTS");

    std::string hashstr;

    if(usecache)
    {
        // Need the name and version of a module for the cache lookup
        auto minfo = MManager().ModuleKeyInfo(key);

        // Create a hash for the lookup
        auto hash = MakeHash(HashType::Hash128,
                             minfo.name,
                             minfo.version,
                             deriv, wfn, bs1, bs2);

        hashstr = hash_to_string(hash);
        out.Debug("Going to lookup one-electron integrals with hash %?\n", hashstr);

        if(Cache().Count(hashstr))
        {
            out.Debug("Integrals were found in the cache. Returning\n");
            return Cache().Get<ReturnType>(hashstr);
        }
    }


    out.Debug("Integrals not found or cache is not being used. Calculating\n");

    const size_t nshell1 = bs1.NShell();
    const size_t nshell2 = bs2.NShell();
    const size_t nfunc1 = bs1.NFunctions();
    const size_t nfunc2 = bs2.NFunctions();

    // actually create the module
    auto mod = CreateChild<OneElectronIntegral>(key);
    mod->Initialize(deriv, wfn, bs1, bs2);
    const unsigned int ncomp = mod->NComponents(); 

    const size_t bufsize = ncomp * bs1.MaxNFunctions() * bs2.MaxNFunctions();
    std::unique_ptr<double []> buffer(new double[bufsize]);

    // vector of ncomp elements, each created with the proper size
    std::vector<MatrixXd> mats(ncomp, MatrixXd(nfunc1, nfunc2));

    for(size_t n1 = 0; n1 < nshell1; n1++)
    {
        const auto & sh1 = bs1.Shell(n1);
        const size_t rowstart = bs1.ShellStart(n1);

        for(size_t n2 = 0; n2 < nshell2; n2++)
        {
            const auto & sh2  = bs2.Shell(n2);
            const size_t colstart = bs2.ShellStart(n2);

            // calculate
            size_t ncalc = mod->Calculate(n1, n2, buffer.get(), bufsize);

            // iterate and fill in the matrix
            AOIterator<2> aoit({sh1, sh2}, false);

            // make sure the right number of integrals was returned
            if(ncalc != aoit.NFunctions() * ncomp)
                throw GeneralException("Bad number of integrals returned",
                                       "ncalc", ncalc, "expected", ncomp * aoit.NFunctions());

            do {
                const size_t i = rowstart+aoit.ShellFunctionIdx<0>();
                const size_t j = colstart+aoit.ShellFunctionIdx<1>();

                for(unsigned int c = 0; c < ncomp; c++)
                    mats[c](i,j) = buffer[aoit.TotalIdx()];

            } while(aoit.Next());
        }
    }


    // convert to the appropriate type
    ReturnType ret;

    for(unsigned int i = 0; i < ncomp; i++)
        ret.push_back(std::make_shared<EigenMatrixImpl>(std::move(mats[i])));

    // Put in cache
    if(usecache)
        Cache().Set(hashstr, ret); // note - ret is a vector of shared_ptr, so this is ok
 
    return ret;
}

