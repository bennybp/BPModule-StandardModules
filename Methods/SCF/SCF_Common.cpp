#include <pulsar/output/OutputStream.hpp>
#include <pulsar/system/AOIterator.hpp>
#include "Methods/SCF/SCF_Common.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace pulsar::modulemanager;
using namespace pulsar::modulebase;
using namespace pulsar::system;
using namespace pulsar::math;



namespace pulsarmethods{


MatrixXd
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

    for(size_t n1 = 0; n1 < nshell; n1++)
    {
        const auto & sh1 = bs.Shell(n1);
        const size_t rowstart = bs.ShellStart(n1);
        const size_t nfunc1 = sh1.NFunctions();

        for(size_t n2 = 0; n2 <= n1; n2++)
        {
            const auto & sh2  = bs.Shell(n2);
            const size_t colstart = bs.ShellStart(n2);
            const size_t nfunc2 = sh2.NFunctions();

            // calculate
            size_t ncalc = mod->Calculate(0, n1, n2, b.data(), maxnfunc2);

            // iterate and fill in the matrix
            AOIterator<2> aoit({sh1, sh2}, false);

            do {
                const size_t i = rowstart+aoit.ShellFunctionIdx<0>();
                const size_t j = colstart+aoit.ShellFunctionIdx<1>();

                mat(i,j) = mat(j, i) = b[aoit.TotalIdx()];
            } while(aoit.Next());
        }
    }

    return mat;
}



std::vector<double>
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


BlockByIrrepSpin<MatrixXd> FormDensity(const BlockByIrrepSpin<MatrixXd> & Cmat,
                                       const BlockByIrrepSpin<VectorXd> & occ)
{
    BlockByIrrepSpin<MatrixXd> Dmat;

    for(auto s : Cmat.GetSpins(Irrep::A))
    {
        const MatrixXd & c = Cmat.Get(Irrep::A, s);
        const VectorXd & o = occ.Get(Irrep::A, s);

        MatrixXd d(c.rows(), c.cols());
        for(size_t i = 0; i < c.rows(); i++)
        for(size_t j = 0; j < c.cols(); j++)
        {
            d(i,j) = 0.0;
            for(size_t m = 0; m < o.size(); m++)
                d(i,j) += o(m) * c(i,m) * c(j,m);
        }
        Dmat.Take(Irrep::A, s, std::move(d));
    }

    return Dmat;
}



BlockByIrrepSpin<Eigen::VectorXd> FindOccupations(size_t nelec)
{
    BlockByIrrepSpin<Eigen::VectorXd> occ;


    if(nelec %2 == 0)
    {
        size_t ndocc = nelec/2;
        VectorXd docc(ndocc);
        for(size_t i = 0; i < ndocc; i++)
            docc(i) = 2.0;
        occ.Take(Irrep::A, 0, std::move(docc));
    }
    else
    {
        size_t nbetaocc = nelec/2; // integer division
        size_t nalphaocc = nelec - nbetaocc;

        VectorXd alphaocc(nalphaocc), betaocc(nbetaocc);
        for(size_t i = 0; i < nalphaocc; i++) alphaocc(i) = 1.0;
        for(size_t i = 0; i < nbetaocc; i++) betaocc(i) = 1.0;

        occ.Take(Irrep::A,  1, std::move(alphaocc));
        occ.Take(Irrep::A, -1, std::move(betaocc));
    }

    return occ;
}



SimpleMatrixD EigenToSimpleMatrix(const Eigen::MatrixXd & m)
{
    // eigen stores in column major by default
    SimpleMatrixD s(m.rows(), m.cols());
    for(size_t i = 0; i < m.rows(); i++)
    for(size_t j = 0; j < m.cols(); j++)
        s(i,j) = m(i,j);
    return s;
}

SimpleVectorD EigenToSimpleVector(const Eigen::VectorXd & v)
{
    return SimpleVectorD(v.size(), v.data());
}

Eigen::MatrixXd SimpleMatrixToEigen(const pulsar::math::SimpleMatrixD & m)
{
    using Eigen::Dynamic;
    using Eigen::RowMajor;

    //! \todo can't use map because of const issues?
    MatrixXd em(m.NRows(), m.NCols());
    for(size_t i = 0; i < m.NRows(); i++)
    for(size_t j = 0; j < m.NCols(); j++)
        em(i,j) = m(i,j);
    return em;
}

Eigen::VectorXd SimpleVectorToEigen(const pulsar::math::SimpleVectorD & v)
{
    //! \todo can't use map because of const issues?
    Eigen::VectorXd ret(v.Size());
    std::copy(v.Data(), v.Data() + v.Size(), ret.data());
    return ret;
}


}//End namespace
