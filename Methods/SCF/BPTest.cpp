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
    


BlockByIrrepSpin<MatrixXd> BuildFock(const BlockByIrrepSpin<MatrixXd> & Dmat,
                                     const MatrixXd & Hcore,
                                     const std::vector<double> & eri)
{
    BlockByIrrepSpin<MatrixXd> Fmat;

    for(auto s : Dmat.GetSpins(Irrep::A))
    {
        MatrixXd F = Hcore;
        auto D = Dmat.Get(Irrep::A, s);

        size_t nao1 = Hcore.rows();
        size_t nao2 = Hcore.cols();

        for(size_t mu = 0; mu < nao1; mu++)
        for(size_t nu = 0; nu < nao2; nu++)
        {
            for(size_t lambda = 0; lambda < nao1; lambda++)
            for(size_t sigma = 0; sigma < nao2; sigma++)
            {
                size_t mnls = INDEX4(mu, nu, lambda, sigma);
                size_t mlns = INDEX4(mu, lambda, nu, sigma);
                F(mu, nu) += 0.5*D(lambda, sigma) * (2*eri.at(mnls)-eri.at(mlns));
            }
        }

        Fmat.Take(Irrep::A, s, std::move(F));
    }

    return Fmat;
}



double CalculateEnergy(const BlockByIrrepSpin<MatrixXd> & Dmat,
                       const BlockByIrrepSpin<MatrixXd> & Fmat,
                       const MatrixXd & Hcore,
                       OutputStream & out)
{
    // calculate the energy 
    double energy = 0.0;
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
    energy = oneelectron + twoelectron;

    out.Output("        One electron: %12.8e\n", oneelectron);
    out.Output("        Two electron: %12.8e\n", twoelectron);
    out.Output("    Total Electronic: %12.8e\n", energy);

    return energy;
}

    
BPTest::DerivReturnType BPTest::Deriv_(size_t order, const Wavefunction & wfn)
{
    if(order != 0)
        throw NotYetImplementedException("Test with deriv != 0");

    if(!wfn.system)
        throw GeneralException("System is not set!");

    // get the basis set
    const System & sys = *(wfn.system);
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
    // These represent the "current" cmatrix, etc
    //////////////////////////////////////////////////////
    BlockByIrrepSpin<MatrixXd> Cmat, Dmat, Fmat;
    BlockByIrrepSpin<VectorXd> epsilon;
    BlockByIrrepSpin<VectorXd> occ;
    double energy = 0;
 
    //////////////////////////
    // Initial Guess
    //////////////////////////
    if(!wfn.cmat) // c-matrix hasn't been set in the passwd wfn
    {
        out.Debug("Don't have C-matrices set. Will call initial guess module\n");
        if(!Options().Has("KEY_INITIAL_GUESS"))
            throw GeneralException("Missing initial guess module when I don't have a C-matrix");
        auto mod_iguess = CreateChildFromOption<EnergyMethod>("KEY_INITIAL_GUESS");
        auto iguess_ret = mod_iguess->Energy(wfn);
        const Wavefunction & guesswfn = iguess_ret.first;

        Cmat = guesswfn.cmat->TransformType<Eigen::MatrixXd>(SimpleMatrixToEigen);
        occ = guesswfn.occupations->TransformType<Eigen::VectorXd>(SimpleVectorToEigen);
        epsilon = guesswfn.epsilon->TransformType<Eigen::VectorXd>(SimpleVectorToEigen);
    }
    else
    {
        out.Debug("Using initial wavefunction as a starting point");
 
        // c-matrices have been set. Make sure we have occupations, etc, as well
        if(!wfn.occupations)
            throw GeneralException("Missing Occupations");
        if(!wfn.epsilon)
            throw GeneralException("Missing Epsilon");

        Cmat = wfn.cmat->TransformType<Eigen::MatrixXd>(SimpleMatrixToEigen);
        occ = wfn.occupations->TransformType<Eigen::VectorXd>(SimpleVectorToEigen);
        epsilon = wfn.epsilon->TransformType<Eigen::VectorXd>(SimpleVectorToEigen);
    }


    ///////////////////////////////////////////
    // Load the one electron integral matrices
    // (and nuclear repulsion)
    ///////////////////////////////////////////
    // Nuclear repulsion
    auto mod_nuc_rep = CreateChildFromOption<SystemIntegral>("KEY_NUC_REPULSION");
    double nucrep;
    mod_nuc_rep->Calculate(0, *wfn.system, &nucrep, 1); //! \todo ignore return value?

    /////////////////////// 
    // Overlap
    auto mod_ao_overlap = CreateChildFromOption<OneElectronIntegral>("KEY_AO_OVERLAP");
    mod_ao_overlap->SetBases(*wfn.system, bstag, bstag);
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


    //////////////////////////// 
    // One-electron hamiltonian
    auto mod_ao_core = CreateChildFromOption<OneElectronIntegral>("KEY_AO_COREBUILD");
    mod_ao_core->SetBases(*wfn.system, bstag, bstag);
    MatrixXd Hcore = FillOneElectronMatrix(mod_ao_core, bs);


    
    //////////////////////////////////////////////////////////////////
    // Actual SCF Procedure
    //////////////////////////////////////////////////////////////////
    // Calculate the initial Density
    Dmat = FormDensity(Cmat, occ);

    // Initial fock matrix is just the core Hamiltonain
    for(const int s : Dmat.GetSpins(Irrep::A))
        Fmat.Set(Irrep::A, s, Hcore);

    // initial energy
    out.Output("Initial guess\n");
    double elec_energy = CalculateEnergy(Dmat, Fmat, Hcore, out);
    energy = elec_energy + nucrep;
    out.Output("   Nuclear Repulsion: %12.8e\n", nucrep);
    out.Output("               Total: %12.8e\n", energy);


    // Storing the results of the previous iterations
    BlockByIrrepSpin<MatrixXd> lastDmat, lastCmat;
    BlockByIrrepSpin<VectorXd> lastocc, lastepsilon;
    double lastenergy;


    // obtain the options
    double etol = Options().Get<double>("E_TOLERANCE");
    size_t maxniter = Options().Get<double>("MAX_ITER");

    size_t iter = 0;
    do
    {
        iter++; 

        auto mod_iter = CreateChildFromOption<EnergyMethod>("KEY_SCF_ITERATOR");

        Wavefunction newwfn;
        newwfn.system = wfn.system;
        newwfn.cmat = std::make_shared<const IrrepSpinMatrixD>(Cmat.TransformType<SimpleMatrixD>(EigenToSimpleMatrix));
        newwfn.occupations = std::make_shared<const IrrepSpinVectorD>(occ.TransformType<SimpleVectorD>(EigenToSimpleVector));
        newwfn.epsilon = std::make_shared<const IrrepSpinVectorD>(occ.TransformType<SimpleVectorD>(EigenToSimpleVector));

        auto energy_ret = mod_iter->Energy(newwfn);

        // store what we did before
        lastenergy = energy;
        //lastDmat = std::move(Dmat);
        lastCmat = std::move(Cmat);
        lastocc = std::move(occ);
        lastepsilon = std::move(epsilon);


        // New energy and density
        energy = energy_ret.second;
        const Wavefunction & retwfn = energy_ret.first;
        Cmat = retwfn.cmat->TransformType<MatrixXd>(SimpleMatrixToEigen);
        occ = retwfn.occupations->TransformType<VectorXd>(SimpleVectorToEigen);
        epsilon = retwfn.epsilon->TransformType<VectorXd>(SimpleVectorToEigen);

        // new density
        //Dmat = FormDensity(Cmat, occ);

        out.Output("Iteration %?\n", iter);
        out.Output("       Total energy: %12.8e\n", energy);

        out.Output(" Difference from last step: %12.8e\n", energy - lastenergy);

    } while(fabs(energy-lastenergy) > etol && iter < maxniter);

    // What are we returning 
    Wavefunction newwfn;
    newwfn.system = wfn.system;
    newwfn.cmat = std::make_shared<const IrrepSpinMatrixD>(Cmat.TransformType<SimpleMatrixD>(EigenToSimpleMatrix));
    newwfn.occupations = std::make_shared<const IrrepSpinVectorD>(occ.TransformType<SimpleVectorD>(EigenToSimpleVector));
    newwfn.epsilon = std::make_shared<const IrrepSpinVectorD>(epsilon.TransformType<SimpleVectorD>(EigenToSimpleVector));

    return {std::move(newwfn), {energy}};
}
    

}//End namespace
