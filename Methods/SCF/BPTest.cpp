#include "Methods/SCF/BPTest.hpp"
#include "Methods/SCF/SCF_Common.hpp"

#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOIterator.hpp>
#include <pulsar/math/Cast.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;

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

                    AOIterator<4> aoit({sh1, sh2, sh3, sh4}, false);

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

    if(!iwfn.GetSystem())
        throw GeneralException("System is not set!");

    //if(!iwfn.cmat)
    //    throw GeneralException("C matrix is not set!");

    // get the basis set
    const System & sys = *(iwfn.GetSystem());
    std::string bstag = Options().Get<std::string>("BASIS_SET");

    out.Output("Obtaining basis set %? from system\n", bstag);
    const BasisSet bs = sys.GetBasisSet(bstag);
    size_t nao = bs.NFunctions();
    size_t nshell = bs.NShell();
    size_t maxnfunc = bs.MaxNFunctions();
    size_t maxnfunc2 = maxnfunc * maxnfunc;

    out.Output("NAO: %? nshell: %?\n", nao, nshell);
    bs.Print(out);
  
    //////////////////////////////////////////////////////
    // Storage of eigen matrices, etc, by irrep and spin 
    //////////////////////////////////////////////////////
    BlockByIrrepSpin<MatrixXd> Cmat, Dmat, Fmat;
    BlockByIrrepSpin<VectorXd> epsilon;
    BlockByIrrepSpin<VectorXd> occ;
    double energy = 0;
 
    //////////////////////////
    // Initial Guess
    //////////////////////////
    if(!InitialWfn().HasCMat()) // c-matrix hasn't been set
    {
        out.Debug("Don't have C-matrices set. Will call initial guess module\n");
        if(!Options().Has("KEY_INITIAL_GUESS"))
            throw GeneralException("Missing initial guess module when I don't have a C-matrix");
        auto mod_iguess = CreateChildFromOption<EnergyMethod>("KEY_INITIAL_GUESS");
        mod_iguess->Energy();

        // load the cmatrix and occupations from there
        const auto & guesswfn = mod_iguess->FinalWfn();

        Cmat = guesswfn.GetCMat()->TransformType<Eigen::MatrixXd>(SimpleMatrixToEigen);
        occ = guesswfn.GetOccupations()->TransformType<Eigen::VectorXd>(SimpleVectorToEigen);
        epsilon = guesswfn.GetEpsilon()->TransformType<Eigen::VectorXd>(SimpleVectorToEigen);
    }
    else
    {
        out.Debug("Using initial wavefunction as a starting point");
 
        // c-matrices have been set. Make sure we have occupations, etc, as well
        if(!InitialWfn().HasOccupations())
            throw GeneralException("Missing Occupations");
        if(!InitialWfn().HasEpsilon())
            throw GeneralException("Missing Epsilon");

        const auto & iwfn = InitialWfn();
        Cmat = iwfn.GetCMat()->TransformType<Eigen::MatrixXd>(SimpleMatrixToEigen);
        occ = iwfn.GetOccupations()->TransformType<Eigen::VectorXd>(SimpleVectorToEigen);
        epsilon = iwfn.GetEpsilon()->TransformType<Eigen::VectorXd>(SimpleMatrixToEigen);
    }


    ///////////////////////////////////////////
    // Load the one electron integral matrices
    // (and nuclear repulsion)
    ///////////////////////////////////////////
    // Nuclear repulsion
    auto mod_nuc_rep = CreateChildFromOption<SystemIntegral>("KEY_NUC_REPULSION");
    double nucrep;
    size_t n = mod_nuc_rep->Calculate(0, &nucrep, 1);

    /////////////////////// 
    // Overlap
    auto mod_ao_overlap = CreateChildFromOption<OneElectronIntegral>("KEY_AO_OVERLAP");
    mod_ao_overlap->SetBases(bstag, bstag);
    MatrixXd overlap_mat = FillOneElectronMatrix(mod_ao_overlap, bs);

    // diagonalize the overlap
    SelfAdjointEigenSolver<MatrixXd> esolve(overlap_mat);
    MatrixXd s_evec = esolve.eigenvectors();
    VectorXd s_eval = esolve.eigenvalues();

    // not sure an easier way to do this
    for(size_t i = 0; i < s_eval.size(); i++)
        s_eval(i) = 1.0/sqrt(s_eval(i));

    // the S^(-1/2) matrix
    MatrixXd S12 = s_evec * s_eval.asDiagonal() * s_evec.transpose();

    /////////////////////// 
    // Nuclear Attraction
    auto mod_ao_nucatt = CreateChildFromOption<OneElectronIntegral>("KEY_AO_NUCATT");
    mod_ao_nucatt->SetBases(bstag, bstag);
    MatrixXd nucatt_mat = FillOneElectronMatrix(mod_ao_nucatt, bs);

    /////////////////////// 
    // Kinetic Energy
    auto mod_ao_kinetic = CreateChildFromOption<OneElectronIntegral>("KEY_AO_KINETIC");
    mod_ao_kinetic->SetBases(bstag, bstag);
    MatrixXd kinetic_mat = FillOneElectronMatrix(mod_ao_kinetic, bs);

    /////////////////////////
    // Load the ERI to core
    /////////////////////////
    auto mod_ao_eri = CreateChildFromOption<TwoElectronIntegral>("KEY_AO_ERI");
    mod_ao_eri->SetBases(bstag, bstag, bstag, bstag);
    std::vector<double> eri = FillTwoElectronVector(mod_ao_eri, bs);


    



    //////////////////////////////////////////////////////////////////
    // Procedure
    //////////////////////////////////////////////////////////////////
    // 1. Form the core
    MatrixXd Hcore = nucatt_mat + kinetic_mat;  // H = T + V

    // 2. Initial fock matrix
    MatrixXd F0 = S12.transpose() * Hcore * S12.transpose();
    SelfAdjointEigenSolver<MatrixXd> fsolve(F0);
    MatrixXd C0 = fsolve.eigenvectors();
    VectorXd e0 = fsolve.eigenvalues();

    // Tranform C0
    C0 = S12*C0;


    // The initial C matrix is the same for all spins
    // Use the irrep/spin from occupations
    for(auto s : occ.GetSpins(Irrep::A))
        Cmat.Set(Irrep::A, s, C0);


    // Calculate the initial Density
    for(auto s : Cmat.GetSpins(Irrep::A))
    {
        const auto & c = Cmat.Get(Irrep::A, s);
        const auto & o = occ.Get(Irrep::A, s);

        MatrixXd d(c.rows(), c.cols());

        for(size_t i = 0; i < c.rows(); i++)
        for(size_t j = 0; j < c.cols(); j++)
        {
            d(i,j) = 0.0;
            for(size_t m = 0; m < o.size(); m++)
                d(i,j) += o(m) * c(i,m) * c(j,m);
        }

        out << "Initial Dmat:\n" << d << "\n";
        Dmat.Take(Irrep::A, s, std::move(d));
    }

    // initial energy
    energy = 0;
    for(auto s : Dmat.GetSpins(Irrep::A))
    {
        const auto & d = Dmat.Get(Irrep::A, s);
        for(size_t i = 0; i < d.rows(); i++)
        for(size_t j = 0; j < d.cols(); j++)
            energy += d(i,j) * Hcore(i,j);
    }


    out.Output("Initial guess energy: %12.8e\n", energy);

    double lastenergy = energy;

    BlockByIrrepSpin<MatrixXd> lastDmat = Dmat;

    size_t iter = 0;
    //do
    for(int i = 0; i < 50; i++)
    {
        iter++; 

        // New fock matrix
        for(auto s : Dmat.GetSpins(Irrep::A))
        {
            MatrixXd F = Hcore;
            auto lastD = lastDmat.Get(Irrep::A, s);

            for(size_t mu = 0; mu < nao; mu++)
            for(size_t nu = 0; nu < nao; nu++)
            {
                for(size_t lambda = 0; lambda < nao; lambda++)
                for(size_t sigma = 0; sigma < nao; sigma++)
                {
                    size_t mnls = INDEX4(mu, nu, lambda, sigma);
                    size_t mlns = INDEX4(mu, lambda, nu, sigma);
                    F(mu, nu) += 0.5*lastD(lambda, sigma) * (2*eri.at(mnls)-eri.at(mlns));
                }
            }


            // store this
            Fmat.Take(Irrep::A, s, std::move(F));
        }
    
        // New density
        for(auto s : Fmat.GetSpins(Irrep::A))
        {
            MatrixXd Fprime = S12.transpose() * Fmat.Get(Irrep::A, s) * S12;

            SelfAdjointEigenSolver<MatrixXd> fsolve(Fprime);
            MatrixXd c = fsolve.eigenvectors();
            VectorXd e = fsolve.eigenvalues();
            c = S12*c;

            const auto & o = occ.Get(Irrep::A, s);

            MatrixXd d(c.rows(), c.cols());
            for(size_t i = 0; i < c.rows(); i++)
            for(size_t j = 0; j < c.cols(); j++)
            {
                d(i,j) = 0.0;
                for(size_t m = 0; m < o.size(); m++)
                    d(i,j) += o(m) * c(i,m) * c(j,m);
            }

            Dmat.Take(Irrep::A, s, std::move(d));
            Cmat.Take(Irrep::A, s, std::move(c));
            epsilon.Take(Irrep::A, s, std::move(e));
        }

        // calculate the energy 
        double oneelectron = 0.0;
        double twoelectron = 0.0;
        for(auto s : Dmat.GetSpins(Irrep::A))
        {
            const auto & d = Dmat.Get(Irrep::A, s);
            const auto & f = Fmat.Get(Irrep::A, s);
            for(size_t i = 0; i < d.rows(); i++)
            for(size_t j = 0; j < d.cols(); j++)
            {
                oneelectron += d(i,j) * Hcore(i,j);
                twoelectron += 0.5 * d(i,j) * f(i,j);
            }

        }

        twoelectron -= 0.5*oneelectron;
        energy = oneelectron + twoelectron + nucrep;

        out.Output("Iteration %? - Energy = %12.8e\n", iter, energy);
        out.Output("        One electron: %12.8e\n", oneelectron);
        out.Output("        Two electron: %12.8e\n", twoelectron);
        out.Output("   Nuclear Repulsion: %12.8e\n", nucrep);
        out.Output("               Total: %12.8e\n", oneelectron + twoelectron + nucrep);

        lastDmat = Dmat;
    }//while(fabs(energy-lastenergy) > 1e-6);

    // create the occupations and other data that will eventually
    // be saved
    //std::shared_ptr<IrrepSpinVectorD> occ(new IrrepSpinVectorD);
    //std::shared_ptr<IrrepSpinVectorD> epsilon(new IrrepSpinVectorD);
    //std::shared_ptr<IrrepSpinMatrixD> cmat(new IrrepSpinMatrixD);

    // store the final information
    //FinalWfn().SetOccupations(irrepspinocc);

    return {0.0};
}
    

}//End namespace
