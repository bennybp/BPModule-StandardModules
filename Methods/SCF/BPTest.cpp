#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/system/AOIterator.hpp>
#include <pulsar/math/Cast.hpp>
#include <eigen3/Eigen/Dense>
#include "Methods/SCF/BPTest.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using Eigen::RowMajor;
using Eigen::ColMajor;
using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::SelfAdjointEigenSolver;

using namespace pulsar::modulemanager;
using namespace pulsar::modulebase;
using namespace pulsar::exception;
using namespace pulsar::system;
using namespace pulsar::math;
using namespace pulsar::output;
using namespace pulsar::datastore;

template class pulsar::system::AOIterator<2>;

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

        for(size_t n2 = 0; n2 <= n1; n2++)
        {
            const auto & sh2  = bs.Shell(n2);
            const size_t nfunc2 = sh2.NFunctions();

            // calculate
            size_t ncalc = mod->Calculate(0, n1, n2, b.data(), maxnfunc2); 

            // iterate and fill in the matrix
            AOIterator<2> aoit({sh1, sh2});

            do { 
                const size_t i = rowstart+aoit.ShellFunctionIdx<0>();
                const size_t j = colstart+aoit.ShellFunctionIdx<1>();

                mat(i,j) = mat(j, i) = b[aoit.TotalIdx()];
            } while(aoit.Next());

            colstart += nfunc2;
        }

        rowstart += nfunc1;
    }

    return mat;
}


#define INDEX2(i,j) (  (j > i) ? (j*(j+1))/2 + i : (i*(i+1))/2 + j )
#define INDEX4(i,j,k,l)  ( INDEX2(k,l) > INDEX2(i,j) ?  (INDEX2(k,l)*(INDEX2(k,l)+1))/2 + INDEX2(i,j) : (INDEX2(i,j)*(INDEX2(i,j)+1))/2 + INDEX2(k,l) )

static std::vector<double>
FillTwoElectronVector(ModulePtr<TwoElectronIntegral> & mod,
                      const BasisSet & bs)
{
    const size_t nao = bs.NFunctions();
    const size_t nshell = bs.NShell();
    const size_t maxnfunc = bs.MaxNFunctions();

    const size_t nao12 = (nao*(nao+1))/2;
    const size_t nao1234 = (nao12*(nao12+1))/2;
    const size_t bufsize = maxnfunc*maxnfunc*maxnfunc*maxnfunc;

    std::vector<double> eri(nao1234);
    std::vector<double> eribuf(bufsize);

    size_t i_start = 0;
    for(size_t i = 0; i < nshell; i++)
    {
        const auto & sh1 = bs.Shell(i);
        const size_t ng1 = sh1.NGeneral();
        size_t j_start = 0;

        for(size_t j = 0; j <= i; j++)
        {
            const auto & sh2 = bs.Shell(j);
            const size_t ng2 = sh2.NGeneral();
            size_t k_start = 0;

            for(size_t k = 0; k < nshell; k++)
            {
                const auto & sh3 = bs.Shell(k);
                const size_t ng3 = sh3.NGeneral();
                size_t l_start = 0;

                for(size_t l = 0; l <= k; l++)
                {
                    if(INDEX2(k,l) > INDEX2(i,j))
                        continue;

                    const auto & sh4 = bs.Shell(l);
                    const size_t ng4 = sh4.NGeneral();

                    uint64_t ncalc = mod->Calculate(0, i, j, k, l, eribuf.data(), bufsize); 

                    AOIterator<4> aoit({sh1, sh2, sh3, sh4});

                    do { 
                        const size_t full_i = i_start+aoit.ShellFunctionIdx<0>();
                        const size_t full_j = j_start+aoit.ShellFunctionIdx<1>();
                        const size_t full_k = k_start+aoit.ShellFunctionIdx<2>();
                        const size_t full_l = l_start+aoit.ShellFunctionIdx<3>();

                        eri.at(INDEX4(full_i, full_j, full_k, full_l)) = eribuf.at(aoit.TotalIdx());
                    } while(aoit.Next());

                    l_start += sh4.NFunctions();   
                }
                k_start += sh3.NFunctions();
            }
            j_start += sh2.NFunctions();
        }
        i_start += sh1.NFunctions();
    }

    return std::move(eri);
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
    
    double nelec_d = sys.GetNElectrons();
    if(!IsInteger(nelec_d))
        throw GeneralException("Can't handle non-integer occupations", "nelectrons", nelec_d);
    size_t nelec = numeric_cast<size_t>(nelec_d);
    out.Output("Number of electrons: %?\n", nelec);



    // Step 1: Nuclear repulsion
    auto mod_nuc_rep = CreateChildFromOption<SystemIntegral>("KEY_NUC_REPULSION");
    double nucrep;
    size_t n = mod_nuc_rep->Calculate(0, &nucrep, 1);
    out.Output("Nuclear repulsion: %12.8e\n", nucrep);


    // Step 2: Overlap
    auto mod_ao_overlap = CreateChildFromOption<OneElectronIntegral>("KEY_AO_OVERLAP");
    mod_ao_overlap->SetBases(bstag, bstag);
    MatrixXd overlap_mat = FillOneElectronMatrix(mod_ao_overlap, bs);

    out << "\nOverlap matrix\n";
    out << overlap_mat << "\n";


    // diagonalize the overlap
    SelfAdjointEigenSolver<MatrixXd> esolve(overlap_mat);
    MatrixXd s_evec = esolve.eigenvectors();
    VectorXd s_eval = esolve.eigenvalues();


    out << "\nEigenvectors of the overlap matrix\n";
    out << s_evec << "\n";

    out << "\nEigenvalues of the overlap matrix\n";
    out << s_eval << "\n";

    // not sure an easier way to do this
    for(size_t i = 0; i < s_eval.size(); i++)
        s_eval(i) = 1.0/sqrt(s_eval(i));


    MatrixXd S12 = s_evec * s_eval.asDiagonal() * s_evec.transpose();

    
    out << "\nS^(1/2) Matrix\n";
    out << S12 << "\n";


    // Step 3: Nuclear Attraction
    auto mod_ao_nucatt = CreateChildFromOption<OneElectronIntegral>("KEY_AO_NUCATT");
    mod_ao_nucatt->SetBases(bstag, bstag);
    MatrixXd nucatt_mat = FillOneElectronMatrix(mod_ao_nucatt, bs);

    out << "\nElectron-Nuclear attraction\n";
    out << nucatt_mat << "\n";


    // Step 4: Kinetic Energy
    auto mod_ao_kinetic = CreateChildFromOption<OneElectronIntegral>("KEY_AO_KINETIC");
    mod_ao_kinetic->SetBases(bstag, bstag);
    MatrixXd kinetic_mat = FillOneElectronMatrix(mod_ao_kinetic, bs);

    out << "\nKinetic energy matrix\n";
    out << kinetic_mat << "\n";


    // Step 5: ERI
    auto mod_ao_eri = CreateChildFromOption<TwoElectronIntegral>("KEY_AO_ERI");
    mod_ao_eri->SetBases(bstag, bstag, bstag, bstag);
    std::vector<double> eri = FillTwoElectronVector(mod_ao_eri, bs);




    //////////////////////////////////////////////////////////////////
    // Procedure
    //////////////////////////////////////////////////////////////////
    // 1. Form the core
    MatrixXd Hcore = nucatt_mat + kinetic_mat;  // H = T + V

    out << "\nCore Hamiltonian\n";
    out << Hcore << "\n";

    // 2. Initial Guess density
    MatrixXd F0 = S12.transpose() * Hcore * S12.transpose();

    out << "\nTransformed Fock Matrix\n";
    out << F0 << "\n";


    SelfAdjointEigenSolver<MatrixXd> fsolve(F0);
    MatrixXd C0 = fsolve.eigenvectors();
    VectorXd e0 = fsolve.eigenvalues();

    // Tranform C0
    C0 = S12*C0;

    out << "Initial C matrix\n";
    out << C0 << "\n";

    out << "Initial orbital energies matrix\n";
    out << e0 << "\n";

    MatrixXd occblock = C0.block(0, 0, nao, nelec/2);
    MatrixXd D0 = occblock * occblock.transpose();

    out << "Initial density matrix\n";
    out << D0 << "\n";


    // 3. Calculate the initial guess energy
    double energy = 0;

    for(size_t i = 0; i < nao; i++)
    for(size_t j = 0; j < nao; j++)
        energy += D0(i,j) * ( 2*Hcore(i,j) );

    out.Output("Initial guess energy: %12.8e\n", energy);

    double lastenergy = 0;

    MatrixXd lastD = D0;

    size_t iter = 0;
    //while(fabs(energy-lastenergy) > 1e-6)
    while(iter < 100)
    {
        iter++; 

        // 4. New fock matrix
        MatrixXd F = Hcore;

        for(size_t mu = 0; mu < nao; mu++)
        for(size_t nu = 0; nu < nao; nu++)
        {
            for(size_t lambda = 0; lambda < nao; lambda++)
            for(size_t sigma = 0; sigma < nao; sigma++)
            {
                size_t mnls = INDEX4(mu, nu, lambda, sigma);
                size_t mlns = INDEX4(mu, lambda, nu, sigma);
                F(mu, nu) += lastD(lambda, sigma) * (2*eri.at(mnls)-eri.at(mlns));
            }
        }
    
        //out.Output("Iteration %? - New F:\n", iter);
        //out << F << "\n";

        // 5. New density
        MatrixXd Fprime = S12.transpose() * F * S12;

        SelfAdjointEigenSolver<MatrixXd> fsolve(Fprime);
        MatrixXd C = fsolve.eigenvectors();
        VectorXd e = fsolve.eigenvalues();

        C = S12*C;

        MatrixXd occblock = C.block(0, 0, nao, nelec/2);
        MatrixXd D = occblock * occblock.transpose();

        energy = 0;
        for(size_t i = 0; i < nao; i++)
        for(size_t j = 0; j < nao; j++)
            energy += D(i,j) * ( Hcore(i,j) + F(i,j));


        out.Output("Iteration %? - Energy = %12.8e\n", iter, energy + nucrep);
        lastD = D;
    }




    return {0.0};
}
    

}//End namespace
