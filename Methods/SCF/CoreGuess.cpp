#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/system/AOIterator.hpp>
#include <pulsar/modulebase/All.hpp>
#include <pulsar/math/Cast.hpp>

#include <Eigen/Dense>
#include "Methods/SCF/CoreGuess.hpp"
#include "Methods/SCF/SCFCommon.hpp"

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


CoreGuess::DerivReturnType CoreGuess::deriv_(size_t order, const Wavefunction & wfn)
{
    if(order != 0)
        throw NotYetImplementedException("CoreGuess with deriv != 0");

    // make sure stuff is set in wavefunction
    if(!wfn.system)
        throw GeneralException("System is not set!");


    // get the basis set
    const System & sys = *(wfn.system);
    std::string bstag = options().get<std::string>("BASIS_SET");

    const BasisSet bs = sys.get_basis_set(bstag);


    ///////////////////////////////////////////
    // Load the one electron integral matrices
    // (and nuclear repulsion)
    ///////////////////////////////////////////
    // Nuclear repulsion
    auto mod_nuc_rep = create_child_from_option<SystemIntegral>("KEY_NUC_REPULSION");
    double nucrep;
    mod_nuc_rep->initialize(0, *wfn.system);
    mod_nuc_rep->calculate(&nucrep, 1);

    /////////////////////////////////////
    // The one-electron integral cacher
    /////////////////////////////////////
    auto mod_ao_cache = create_child_from_option<OneElectronMatrix>("KEY_ONEEL_MAT");

    /////////////////////// 
    // Overlap
    /////////////////////// 
    const std::string ao_overlap_key = options().get<std::string>("KEY_AO_OVERLAP");
    auto overlapimpl = mod_ao_cache->calculate(ao_overlap_key, 0, wfn, bs, bs);
    std::shared_ptr<const MatrixXd> overlap_mat = convert_to_eigen(overlapimpl.at(0));  // .at(0) = first (and only) component
    MatrixXd S12 = FormS12(*overlap_mat);

    //////////////////////////// 
    // One-electron hamiltonian
    //////////////////////////// 
    const std::string ao_build_key = options().get<std::string>("KEY_AO_COREBUILD");
    auto Hcoreimpl = mod_ao_cache->calculate(ao_build_key, 0, wfn, bs, bs);
    std::shared_ptr<const MatrixXd> Hcore = convert_to_eigen(Hcoreimpl.at(0));  // .at(0) = first (and only) component


    //////////////////////////
    // Occupations
    //////////////////////////
    double nelec_d = sys.get_n_electrons();
    if(!is_integer(nelec_d))
        throw GeneralException("Can't handle non-integer occupations", "nelectrons", nelec_d);
    size_t nelec = numeric_cast<size_t>(nelec_d);


    // Block some eigen matrices, etc, by irrep and spin
    IrrepSpinMatrixD cmat, dmat;
    IrrepSpinVectorD epsilon, occ;

    // Fill in the occupations
    occ = FindOccupations(nelec);


    // 2. Initial fock matrix
    MatrixXd F0 = S12.transpose() * (*Hcore) * S12.transpose();
    SelfAdjointEigenSolver<MatrixXd> fsolve(F0);
    MatrixXd C0 = fsolve.eigenvectors();
    VectorXd e0 = fsolve.eigenvalues();

    // Tranform C0
    C0 = S12*C0;


    // The initial C matrix is the same for all spins
    // Use the irrep/spin from occupations
    for(auto s : occ.get_spins(Irrep::A))
    {
        cmat.set(Irrep::A, s, std::make_shared<EigenMatrixImpl>(C0));
        epsilon.set(Irrep::A, s, std::make_shared<EigenVectorImpl>(e0));
    }


    // Calculate the initial Density
    for(auto ir : cmat.get_irreps())
    for(auto s : cmat.get_spins(ir))
    {
        std::shared_ptr<const MatrixXd> cptr = convert_to_eigen(cmat.get(ir, s));
        std::shared_ptr<const VectorXd> optr = convert_to_eigen(occ.get(ir, s));
        const auto & o = *optr;
        const auto & c = *cptr;

        MatrixXd d(c.rows(), c.cols());

        for(long i = 0; i < c.rows(); i++)
        for(long j = 0; j < c.cols(); j++)
        {
            d(i,j) = 0.0;
            for(long m = 0; m < o.size(); m++)
                d(i,j) += o(m) * c(i,m) * c(j,m);
        }

        dmat.set(ir, s, std::make_shared<EigenMatrixImpl>(std::move(d)));
    }

    // initial energy
    double energy = 0.0;
    for(auto ir : dmat.get_irreps())
    for(auto s : dmat.get_spins(ir))
    {
        const auto & d = dmat.get(ir, s);
        for(size_t i = 0; i < d->size(0); i++)
        for(size_t j = 0; j < d->size(1); j++)
            energy += d->get_value({i,j}) * (*Hcore)(i,j);
    }


    out.output("Formed initial guess.\n");
    out.output("    Electronic energy:  %16.8e\n", energy);
    out.output("    Nuclear repulsion:  %16.8e\n", nucrep);

    energy += nucrep;
    out.output("         Total energy:  %16.8e\n", energy);

    Wavefunction newwfn;
    newwfn.system = wfn.system;
    newwfn.cmat = std::make_shared<const IrrepSpinMatrixD>(std::move(cmat));
    newwfn.opdm = std::make_shared<const IrrepSpinMatrixD>(std::move(dmat));
    newwfn.occupations = std::make_shared<const IrrepSpinVectorD>(std::move(occ));
    newwfn.epsilon = std::make_shared<const IrrepSpinVectorD>(std::move(epsilon));

    return {std::move(newwfn), {energy}};
}
    

}//End namespace
