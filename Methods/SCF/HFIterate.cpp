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


void HFIterate::Initialize_(const Wavefunction & wfn)
{
    if(!wfn.system)
        throw GeneralException("System is not set!");

    // get the basis set
    const System & sys = *(wfn.system);
    std::string bstag = Options().Get<std::string>("BASIS_SET");

    const BasisSet bs = sys.GetBasisSet(bstag);


    ///////////////////////
    // Overlap
    ///////////////////////
    auto mod_ao_overlap = CreateChildFromOption<OneElectronIntegral>("KEY_AO_OVERLAP");
    mod_ao_overlap->SetBases(sys, bstag, bstag);
    MatrixXd overlap_mat = FillOneElectronMatrix(mod_ao_overlap, bs);

    // diagonalize the overlap and form S^(-1/2)
    SelfAdjointEigenSolver<MatrixXd> esolve(overlap_mat);
    MatrixXd s_evec = esolve.eigenvectors();
    VectorXd s_eval = esolve.eigenvalues();

    // not sure an easier way to do this
    for(size_t i = 0; i < s_eval.size(); i++)
        s_eval(i) = 1.0/sqrt(s_eval(i));

    // the S^(-1/2) matrix
    S12_ = s_evec * s_eval.asDiagonal() * s_evec.transpose();

    initialized_ = true;
}


Wavefunction HFIterate::Next_(const Wavefunction & wfn, const IrrepSpinMatrixD & fmat)
{
    if(!initialized_)
        Initialize_(wfn);

    // The density and C matrix we are returning
    IrrepSpinMatrixD Dmat, Cmat;
    IrrepSpinVectorD epsilon;

    // Diagonalize, etc
    for(auto ir : fmat.GetIrreps())
    for(auto s : fmat.GetSpins(ir))
    {
        MappedConstMatrix mapped_f = MapConstSimpleMatrix(fmat.Get(ir, s));

        MatrixXd Fprime = S12_.transpose() * mapped_f * S12_;

        SelfAdjointEigenSolver<MatrixXd> fsolve(Fprime);
        MatrixXd c = fsolve.eigenvectors();
        VectorXd e = fsolve.eigenvalues();
        c = S12_*c;

        // form the density matrix, and store both in the wfn
        SimpleMatrixD simple_c = EigenToSimpleMatrix(c);
        SimpleMatrixD simple_d = FormDensity(simple_c, wfn.occupations->Get(ir, s));
        SimpleVectorD simple_e = EigenToSimpleVector(e);
        Dmat.Take(ir, s, std::move(simple_d));
        Cmat.Take(ir, s, std::move(simple_c));
        epsilon.Take(ir, s, std::move(simple_e));
    }

    // build the new wavefunction
    Wavefunction newwfn;
    newwfn.system = wfn.system;
    newwfn.cmat = std::make_shared<const IrrepSpinMatrixD>(std::move(Cmat));
    newwfn.opdm = std::make_shared<const IrrepSpinMatrixD>(std::move(Dmat));
    newwfn.occupations = wfn.occupations; // didn't change
    newwfn.epsilon = std::make_shared<const IrrepSpinVectorD>(std::move(epsilon));

    return newwfn;

}


} // close namespace pulsarmethods
