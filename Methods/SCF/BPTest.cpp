#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <eigen3/Eigen/Dense>
#include "Methods/SCF/BPTest.hpp"

using Eigen::MatrixXd;
using Eigen::Map;

using namespace pulsar::modulemanager;
using namespace pulsar::modulebase;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::math;
using namespace pulsar::output;
using namespace pulsar::datastore;

namespace pulsarmethods{



BPTest::BPTest(ID_t id)
    : EnergyMethod(id) { }
    

static MatrixXd
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

    size_t rowstart = 0;
    for(size_t n1 = 0; n1 < nshell; n1++)
    {
        const auto & sh1 = bs.Shell(n1);
        const size_t nfunc1 = sh1.NFunctions();
        size_t colstart = 0;

        for(size_t n2 = 0; n2 < nshell; n2++)
        {
            const auto & sh2  = bs.Shell(n2);
            const size_t nfunc2 = sh2.NFunctions();

            // calculate
            size_t ncalc = mod->Calculate(0, n1, n2, b.data(), maxnfunc2); 

            // go over the general contractions
            size_t growstart = rowstart;
            size_t offset = 0;
            for(size_t g1 = 0; g1 < bs.Shell(n1).NGeneral(); g1++)
            {
                const size_t gnfunc1 = sh1.GeneralNFunctions(g1);
                size_t gcolstart = colstart;
                for(size_t g2 = 0; g2 < sh2.NGeneral(); g2++)
                {
                    const size_t gnfunc2 = sh2.GeneralNFunctions(g2);
                
                    // copy by blocking the matrix and mapping the raw buffer to a matrix
                    mat.block(growstart, gcolstart, gnfunc1, gnfunc2)
                          = Map<MatrixXd>(b.data() + offset, gnfunc1, gnfunc2);

                    offset += gnfunc1*gnfunc2;
                    gcolstart += gnfunc2;
                }
                growstart += gnfunc1;
            }
                
            colstart += nfunc2;
        }

        rowstart += nfunc1;
    }

    return mat;
}

    
std::vector<double> BPTest::Deriv_(size_t order)
{
    if(order != 0)
        throw NotYetImplementedException("Test with deriv != 0");

    // make sure stuff is set in wavefunction
    const Wavefunction & iwfn = InitialWfn();

    if(!iwfn.system)
        throw GeneralException("System is not set!");

    //if(!iwfn.cmat)
    //    throw GeneralException("C matrix is not set!");

    // get the basis set
    const System & sys = *(iwfn.system);
    std::string bstag = Options().Get<std::string>("BASIS_SET");

    out.Output("Obtaining basis set %? from system\n", bstag);
    const BasisSet bs = sys.GetBasisSet(bstag);
    size_t nao = bs.NFunctions();
    size_t nshell = bs.NShell();
    size_t maxnfunc = bs.MaxNFunctions();
    size_t maxnfunc2 = maxnfunc * maxnfunc;

    out.Output("NAO: %? nshell: %?\n", nao, nshell);
    bs.Print(out);
    


    // Step 1: Nuclear repulsion
    auto mod_nuc_rep = CreateChildModuleFromOption<SystemIntegral>("KEY_NUC_REPULSION");
    double nucrep;
    size_t n = mod_nuc_rep->Calculate(0, &nucrep, 1);
    out.Output("Nuclear repulsion: %12.8e\n", nucrep);


    // Step 2: Overlap
    auto mod_ao_overlap = CreateChildModuleFromOption<OneElectronIntegral>("KEY_AO_OVERLAP");
    mod_ao_overlap->SetBases(bstag, bstag);
    MatrixXd overlap_mat = FillOneElectronMatrix(mod_ao_overlap, bs);

    out << "\nOverlap matrix\n";
    out << overlap_mat << "\n";


    // Step 3: Nuclear Attraction
    auto mod_ao_nucatt = CreateChildModuleFromOption<OneElectronIntegral>("KEY_AO_NUCATT");
    mod_ao_nucatt->SetBases(bstag, bstag);
    MatrixXd nucatt_mat = FillOneElectronMatrix(mod_ao_nucatt, bs);

    out << "\nElectron-Nuclear attraction\n";
    out << nucatt_mat << "\n";






    return {0.0};
}
    

}//End namespace
