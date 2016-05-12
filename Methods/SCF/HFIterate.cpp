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


void HFIterate::Initialize_(const System & sys, const std::string & bstag)
{
    const BasisSet bs = sys.GetBasisSet(bstag);

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

    initialized_ = true;
}


Wavefunction HFIterate::Next_(const Wavefunction & wfn, const IrrepSpinMatrixD & fmat)
{
    if(!wfn.system)
        throw GeneralException("System is not set!");

    // get the basis set
    const System & sys = *(wfn.system);
    std::string bstag = Options().Get<std::string>("BASIS_SET");

    if(!initialized_)
        Initialize_(*wfn.system, bstag);

    IrrepSpinMatrixD Cmat;
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

        // convert back to simple matrices
        Cmat.Take(ir, s, EigenToSimpleMatrix(c)); 
        epsilon.Take(ir, s, EigenToSimpleVector(e));
    }

    // set the final wavefunction stuff
    Wavefunction newwfn;
    newwfn.system = wfn.system;
    newwfn.cmat = std::make_shared<const IrrepSpinMatrixD>(std::move(Cmat));
    newwfn.occupations = wfn.occupations; // didn't change
    newwfn.epsilon = std::make_shared<const IrrepSpinVectorD>(std::move(epsilon));

    return newwfn;
}


} // close namespace pulsarmethods
