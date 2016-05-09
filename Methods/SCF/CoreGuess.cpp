#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/system/AOIterator.hpp>
#include <pulsar/math/Cast.hpp>
#include <eigen3/Eigen/Dense>
#include "Methods/SCF/CoreGuess.hpp"
#include "Methods/SCF/SCF_Common.hpp"

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


std::vector<double> CoreGuess::Deriv_(size_t order)
{
    if(order != 0)
        throw NotYetImplementedException("Core guess with deriv != 0");

    // make sure stuff is set in wavefunction
    const Wavefunction & iwfn = InitialWfn();

    if(!iwfn.GetSystem())
        throw GeneralException("System is not set!");
    if(!iwfn.GetOccupations())
        throw GeneralException("Occupations are not set!");

    // get the basis set
    const System & sys = *(iwfn.GetSystem());
    std::string bstag = Options().Get<std::string>("BASIS_SET");

    out.Output("Obtaining basis set %? from system\n", bstag);
    const BasisSet bs = sys.GetBasisSet(bstag);
    size_t nao = bs.NFunctions();
    size_t nshell = bs.NShell();
    size_t maxnfunc = bs.MaxNFunctions();

    double nelec_d = sys.GetNElectrons();
    if(!IsInteger(nelec_d))
        throw GeneralException("Can't handle non-integer occupations", "nelectrons", nelec_d);
    size_t nelec = numeric_cast<size_t>(nelec_d);
    out.Output("Number of electrons: %?\n", nelec);


    // Calculate the Nuclear repulsion
    auto mod_nuc_rep = CreateChildFromOption<SystemIntegral>("KEY_NUC_REPULSION");
    double nucrep;
    size_t n = mod_nuc_rep->Calculate(0, &nucrep, 1);

    // Read an overlap matrix
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


    MatrixXd S12 = s_evec * s_eval.asDiagonal() * s_evec.transpose();

   
    // start the core hamiltonian
    MatrixXd coreH;
 
    // Read the nucelar attraction 
    auto mod_ao_nucatt = CreateChildFromOption<OneElectronIntegral>("KEY_AO_NUCATT");
    mod_ao_nucatt->SetBases(bstag, bstag);
    coreH = FillOneElectronMatrix(mod_ao_nucatt, bs);

    // Kinetic Energy
    auto mod_ao_kinetic = CreateChildFromOption<OneElectronIntegral>("KEY_AO_KINETIC");
    mod_ao_kinetic->SetBases(bstag, bstag);
    coreH += FillOneElectronMatrix(mod_ao_kinetic, bs);

    // transform
    MatrixXd F0 = S12.transpose() * coreH * S12.transpose();

    // diagonalize
    SelfAdjointEigenSolver<MatrixXd> fsolve(F0);
    MatrixXd C0 = fsolve.eigenvectors();
    VectorXd e0 = fsolve.eigenvalues();

    // Transform C0
    C0 = S12*C0;

    // Form the density
    const auto & occ = *(InitialWfn().GetOccupations());
    IrrepSpinMatrixD isdens;

    // calculate the energy as we go
    double energy = nucrep;

    for(const auto & it : occ)
    {
        SimpleMatrixD dens(nao, nao); 
        const auto & o = it.second;

        size_t nocc = it.second.Size();
        for(size_t i = 0; i < nao; i++)
        for(size_t j = 0; j < nao; j++)
        {
            dens(i,j) = 0;
            for(size_t m = 0; m < nocc; m++)
                dens(i,j) += o(m)*C0(i,m)*C0(j,m); 

            energy += dens(i,j) * ( 2*coreH(i,j) );
        }

        isdens.Take(it.first.first, it.first.second, std::move(dens));
    }

    FinalWfn().SetCMat(C0);

    return {energy};
}
    

}//End namespace
