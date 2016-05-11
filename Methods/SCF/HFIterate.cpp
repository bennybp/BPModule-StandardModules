#include "Methods/SCF/HFIterate.hpp"
#include "pulsar/modulebase/All.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;

using namespace pulsar::datastore;
using namespace pulsar::system;
using namespace pulsar::math;
using namespace pulsar::exception;
using namespace pulsar::output;
using namespace pulsar::modulebase;


namespace pulsarmethods {


static
BlockedEigenMatrix BuildFock(const IrrepSpinMatrixD & Dmat,
                             const MatrixXd & Hcore,
                             const std::vector<double> & eri)
{
    BlockedEigenMatrix Fmat;

    for(auto ir : Dmat.GetIrreps())
    for(auto s : Dmat.GetSpins(ir))
    {
        MappedConstMatrix D = MapConstSimpleMatrix(Dmat.Get(ir,s));

        MatrixXd F = Hcore;

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

        Fmat.Take(ir, s, std::move(F));
    }

    return Fmat;
}


double HFIterate::CalculateEnergy(const IrrepSpinMatrixD & Dmat,
                                  const BlockedEigenMatrix & Fmat)
{
    // calculate the energy
    double energy = 0.0;
    double oneelectron = 0.0;
    double twoelectron = 0.0;

    for(auto ir : Dmat.GetIrreps())
    for(auto s : Dmat.GetSpins(ir))
    {
        MappedConstMatrix d = MapConstSimpleMatrix(Dmat.Get(ir, s));
        const auto & f = Fmat.Get(ir, s);
        for(size_t i = 0; i < d.rows(); i++)
        for(size_t j = 0; j < d.cols(); j++)
        {
            oneelectron += d(i,j) * Hcore_(i,j);
            twoelectron += 0.5 * d(i,j) * f(i,j);
        }
    }

    twoelectron -= 0.5*oneelectron;
    energy = oneelectron + twoelectron;

    out.Output("            One electron: %16.8e\n", oneelectron);
    out.Output("            Two electron: %16.8e\n", twoelectron);
    out.Output("        Total Electronic: %16.8e\n", energy);
    out.Output("       Nuclear Repulsion: %16.8e\n", nucrep_);

    energy += nucrep_;
    out.Output("            Total energy: %16.8e\n", energy);

    return energy;
}


void HFIterate::Initialize_(const System & sys, const std::string & bstag)
{
    out.Output("Obtaining basis set %? from system\n", bstag);
    const BasisSet bs = sys.GetBasisSet(bstag);

    /////////////////////////
    // Load the ERI to core
    /////////////////////////
    auto mod_ao_eri = CreateChildFromOption<TwoElectronIntegral>("KEY_AO_ERI");
    mod_ao_eri->SetBases(sys, bstag, bstag, bstag, bstag);
    eri_ = FillTwoElectronVector(mod_ao_eri, bs);


    ///////////////////////////////////////////
    // Load the one electron integral matrices
    // (and nuclear repulsion)
    ///////////////////////////////////////////
    // Nuclear repulsion
    auto mod_nuc_rep = CreateChildFromOption<SystemIntegral>("KEY_NUC_REPULSION");
    size_t n = mod_nuc_rep->Calculate(0, sys, &nucrep_, 1);

    ///////////////////////
    // Overlap
    ///////////////////////
    auto mod_ao_overlap = CreateChildFromOption<OneElectronIntegral>("KEY_AO_OVERLAP");
    mod_ao_overlap->SetBases(sys, bstag, bstag);
    MatrixXd overlap_mat = FillOneElectronMatrix(mod_ao_overlap, bs);

    // diagonalize the overlap
    SelfAdjointEigenSolver<MatrixXd> esolve(overlap_mat);
    MatrixXd s_evec = esolve.eigenvectors();
    VectorXd s_eval = esolve.eigenvalues();

    // not sure an easier way to do this
    for(size_t i = 0; i < s_eval.size(); i++)
        s_eval(i) = 1.0/sqrt(s_eval(i));

    // the S^(-1/2) matrix
    S12_ = s_evec * s_eval.asDiagonal() * s_evec.transpose();


    ////////////////////////////
    // One-electron hamiltonian
    auto mod_ao_core = CreateChildFromOption<OneElectronIntegral>("KEY_AO_COREBUILD");
    mod_ao_core->SetBases(sys, bstag, bstag);
    Hcore_ = FillOneElectronMatrix(mod_ao_core, bs);


    initialized_ = true;
}


HFIterate::DerivReturnType HFIterate::Deriv_(size_t order, const Wavefunction & wfn)
{
    if(order != 0)
        throw NotYetImplementedException("Test with deriv != 0");

    if(!wfn.system)
        throw GeneralException("System is not set!");

    // get the basis set
    const System & sys = *(wfn.system);
    std::string bstag = Options().Get<std::string>("BASIS_SET");


    if(!initialized_)
        Initialize_(*wfn.system, bstag);

    // c-matrices have been set. Make sure we have occupations, etc, as well
    if(!wfn.cmat)
        throw GeneralException("Missing C matrix");
    if(!wfn.occupations)
        throw GeneralException("Missing Occupations");
    if(!wfn.epsilon)
        throw GeneralException("Missing Epsilon");


    //////////////////////////////////////////////////////////////////
    // Procedure
    //////////////////////////////////////////////////////////////////
    // Calculate the initial Density
    IrrepSpinMatrixD Dmat = FormDensity(*wfn.cmat, *wfn.occupations);

    IrrepSpinMatrixD Cmat;
    IrrepSpinVectorD epsilon;

    // Build the fock matrix
    const BlockedEigenMatrix Fmat = BuildFock(Dmat, Hcore_, eri_);

    // Diagonalize, etc
    for(auto ir : Fmat.GetIrreps())
    for(auto s : Fmat.GetSpins(ir))
    {
        MatrixXd Fprime = S12_.transpose() * Fmat.Get(ir, s) * S12_;

        SelfAdjointEigenSolver<MatrixXd> fsolve(Fprime);
        MatrixXd c = fsolve.eigenvectors();
        VectorXd e = fsolve.eigenvalues();
        c = S12_*c;

        // convert back to simple matrices
        Cmat.Take(ir, s, EigenToSimpleMatrix(c)); 
        epsilon.Take(ir, s, EigenToSimpleVector(e));
    }

    // energy
    Dmat = FormDensity(Cmat, *wfn.occupations);
    double energy = CalculateEnergy(Dmat, Fmat);

    // set the final wavefunction stuff
    Wavefunction newwfn;
    newwfn.system = wfn.system;
    newwfn.cmat = std::make_shared<const IrrepSpinMatrixD>(std::move(Cmat));
    newwfn.occupations = wfn.occupations; // didn't change
    newwfn.epsilon = std::make_shared<const IrrepSpinVectorD>(std::move(epsilon));

    return {std::move(newwfn), {energy}};
}


} // close namespace pulsarmethods
