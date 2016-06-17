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
    mod_ao_overlap->Initialize(0, wfn, bs, bs);
    MatrixXd overlap_mat = FillOneElectronMatrix(mod_ao_overlap, bs);

    // diagonalize the overlap and form S^(-1/2)
    SelfAdjointEigenSolver<MatrixXd> esolve(overlap_mat);
    MatrixXd s_evec = esolve.eigenvectors();
    VectorXd s_eval = esolve.eigenvalues();

    // not sure an easier way to do this
    for(int i = 0; i < s_eval.size(); i++)
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
    IrrepSpinMatrixD Cmat;
    IrrepSpinVectorD epsilon;

    // Diagonalize, etc
    for(auto ir : fmat.GetIrreps())
    for(auto s : fmat.GetSpins(ir))
    {
        std::shared_ptr<const MatrixXd> fptr = convert_to_eigen(fmat.Get(ir, s));
        const MatrixXd & f = *fptr;

        MatrixXd Fprime = S12_.transpose() * f * S12_;

        SelfAdjointEigenSolver<MatrixXd> fsolve(Fprime);
        MatrixXd c = fsolve.eigenvectors();
        VectorXd e = fsolve.eigenvalues();
        c = S12_*c;

        Cmat.Take(ir, s, std::make_shared<EigenMatrixImpl>(std::move(c)));
        epsilon.Take(ir, s, std::make_shared<EigenVectorImpl>(std::move(e)));
    }

    // build the density
    IrrepSpinMatrixD Dmat = FormDensity(Cmat, *wfn.occupations);

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
