#include "Common/BasisSetCommon.hpp"

using namespace pulsar::system;
using namespace pulsar::datastore;
using namespace pulsar::output;


// coordinates not used
static BasisShellInfo NormalizeShell_(const BasisShellInfo & shell, const CoordType &)
{
    // a common factor found in normalization
    // (I forgot exactly what it is... :( )
    static const double norm_fac[25] =
    {
    /* l =    0 */  5.56832799683170785             ,
    /* l =    1 */  2.78416399841585392             ,
    /* l =    2 */  4.17624599762378088             ,
    /* l =    3 */  10.4406149940594522             ,
    /* l =    4 */  36.5421524792080827             ,
    /* l =    5 */  164.439686156436372             ,
    /* l =    6 */  904.418273860400048             ,
    /* l =    7 */  5878.71878009260031             ,
    /* l =    8 */  44090.3908506945023             ,
    /* l =    9 */  374768.32223090327              ,
    /* l =   10 */  3560299.06119358106             ,
    /* l =   11 */  37383140.1425326012             ,
    /* l =   12 */  429906111.639124913             ,
    /* l =   13 */  5373826395.48906142             ,
    /* l =   14 */  72546656339.1023291             ,
    /* l =   15 */  1051926516916.98377             ,
    /* l =   16 */  16304861012213.2485             ,
    /* l =   17 */  269030206701518.6               ,
    /* l =   18 */  4708028617276575.5              ,
    /* l =   19 */  87098529419616646.7             ,
    /* l =   20 */  1.69842132368252461e+18         ,
    /* l =   21 */  3.48176371354917545e+19         ,
    /* l =   22 */  7.48579198413072722e+20         ,
    /* l =   23 */  1.68430319642941362e+22         ,
    /* l =   24 */  3.95811251160912202e+23         ,
    };



    // common to all general contractions
    const size_t nprim = shell.n_primitives();
    const double * const alpha = shell.alpha_ptr();

    BasisShellInfo newshell(shell);

    for(size_t n = 0; n < shell.n_general_contractions(); n++)
    {
        const int iam = shell.general_am(n);
        const double am = static_cast<double>(iam);
        const double m = am + 1.5;
        const double m2 = 0.5 * m;

        std::vector<double> coefs = shell.get_coefs(n);

        double sum = 0.0;

        for(size_t i = 0; i < nprim; ++i)
        {
            const double a1 = alpha[i];
            const double c1 = coefs[i];

            for(size_t j = 0; j < nprim; j++)
            {
                const double a2 = alpha[j];
                const double c2 = coefs[j];
                sum += ( c1 * c2 *  pow(a1*a2, m2) ) / ( pow(a1+a2, m) );
            }
        }

        const double norm = 1.0 / sqrt(sum * norm_fac[iam]);

        // apply the rest of the normalization and store
        for (size_t i = 0; i < nprim; ++i)
            newshell.set_coef(n, i, coefs[i] * norm * pow(alpha[i], m2));
    }

    return newshell;
}


std::shared_ptr<BasisSet> NormalizeBasis(CacheData & cache,
                                         OutputStream & out,
                                         const BasisSet & bs)
{
    using bphash::hash_to_string;

    std::string cachekey = std::string("bs:") + hash_to_string(bs.my_hash());
    if(cache.count(cachekey)) // options are unimportant
    {
        auto ret = cache.get<std::shared_ptr<BasisSet>>(cachekey);

        out.debug("Found normalized basis in cache: %? -> %?\n",
                  hash_to_string(bs.my_hash()), hash_to_string(ret->my_hash()));

        return ret;
    }

    auto newbs = std::make_shared<BasisSet>(bs.transform(NormalizeShell_));

    // add to the cache
    cache.set(cachekey, newbs);

    return newbs;
}

