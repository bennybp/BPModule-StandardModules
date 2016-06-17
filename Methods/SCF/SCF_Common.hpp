#ifndef PULSAR_GUARD_SCF__SCF_COMMON_HPP_
#define PULSAR_GUARD_SCF__SCF_COMMON_HPP_

#include <pulsar/modulebase/All.hpp>
#include <pulsar/modulemanager/ModulePtr.hpp>
#include <pulsar/system/BasisSet.hpp>
#include <pulsar/math/BlockByIrrepSpin.hpp>
#include <pulsar/datastore/Wavefunction.hpp>
#include <pulsar/math/TensorImpl.hpp>

#include <Eigen/Dense>

#include <vector>

#define INDEX2(i,j) (  (j > i) ? (j*(j+1))/2 + i : (i*(i+1))/2 + j )
#define INDEX4(i,j,k,l)  ( INDEX2(k,l) > INDEX2(i,j) ?  (INDEX2(k,l)*(INDEX2(k,l)+1))/2 + INDEX2(i,j) : (INDEX2(i,j)*(INDEX2(i,j)+1))/2 + INDEX2(k,l) )

namespace pulsarmethods {


class EigenMatrixImpl : public pulsar::math::TensorImpl<2, double>
{
    public:
        EigenMatrixImpl(const std::shared_ptr<Eigen::MatrixXd> & mat)
            : mat_(mat) { }

        EigenMatrixImpl(const Eigen::MatrixXd & mat)
            : mat_(std::make_shared<Eigen::MatrixXd>(mat)) { }

        EigenMatrixImpl(Eigen::MatrixXd && mat)
            : mat_(std::make_shared<Eigen::MatrixXd>(std::move(mat))) { }

        virtual std::array<size_t, 2> sizes(void) const
        {
            return {static_cast<size_t>((*mat_).rows()),
                    static_cast<size_t>((*mat_).cols())};
        }

        virtual double get_value(std::array<size_t, 2> idx) const
        {
            return (*mat_)(idx[0], idx[1]);
        }

        virtual void set_value(std::array<size_t, 2> idx, double val)
        {
            (*mat_)(idx[0], idx[1]) = val;
        }

        std::shared_ptr<const Eigen::MatrixXd>
        get_matrix(void) const
        {
            return mat_;
        }

    private:
        std::shared_ptr<Eigen::MatrixXd> mat_;
};


class EigenVectorImpl : public pulsar::math::TensorImpl<1, double>
{
    public:
        EigenVectorImpl(const std::shared_ptr<Eigen::VectorXd> & mat)
            : mat_(mat) { }

        EigenVectorImpl(const Eigen::VectorXd & mat)
            : mat_(std::make_shared<Eigen::VectorXd>(mat)) { }

        EigenVectorImpl(Eigen::VectorXd && mat)
            : mat_(std::make_shared<Eigen::VectorXd>(std::move(mat))) { }

        virtual std::array<size_t, 1> sizes(void) const
        {
            return {static_cast<size_t>((*mat_).size())};
        }

        virtual double get_value(std::array<size_t, 1> idx) const
        {
            return (*mat_)(idx[0]);
        }

        virtual void set_value(std::array<size_t, 1> idx, double val)
        {
            (*mat_)(idx[0]) = val;
        }

        std::shared_ptr<const Eigen::VectorXd>
        get_matrix(void) const
        {
            return mat_;
        }

    private:
        std::shared_ptr<Eigen::VectorXd> mat_;
};


inline
std::shared_ptr<const Eigen::MatrixXd>
convert_to_eigen(const std::shared_ptr<const pulsar::math::TensorImpl<2, double>> & ten)
{
    // does the TensorImpl contain an eigen matrix?
    auto test = std::dynamic_pointer_cast<const EigenMatrixImpl>(ten);
    if(test)
        return test->get_matrix();

    // otherwise, convert elementwise
    auto ret = std::make_shared<Eigen::MatrixXd>();

    auto sizes = ten->sizes();
    for(size_t i = 0; i < sizes[0]; i++)
    for(size_t j = 0; j < sizes[0]; j++)
        (*ret)(i,j) = ten->get_value({i,j});
    return ret;
}

inline
std::shared_ptr<const Eigen::VectorXd>
convert_to_eigen(const std::shared_ptr<const pulsar::math::TensorImpl<1, double>> & ten)
{
    // does the TensorImpl contain an eigen matrix?
    auto test = std::dynamic_pointer_cast<const EigenVectorImpl>(ten);
    if(test)
    {
        std::cout << "Here: matrix is an eigen vector\n";
        return test->get_matrix();
    }

    // otherwise, convert elementwize
    auto ret = std::make_shared<Eigen::VectorXd>();

    auto sizes = ten->sizes();
    for(size_t i = 0; i < sizes[0]; i++)
        (*ret)(i) = ten->get_value({i});
    return ret;
}




typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> MappedMatrix;
typedef Eigen::Map<Eigen::VectorXd> MappedVector;
typedef Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> MappedConstMatrix;
typedef Eigen::Map<const Eigen::VectorXd> MappedConstVector;

typedef pulsar::math::BlockByIrrepSpin<Eigen::MatrixXd> BlockedEigenMatrix;
typedef pulsar::math::BlockByIrrepSpin<Eigen::VectorXd> BlockedEigenVector;


Eigen::MatrixXd
FillOneElectronMatrix(pulsar::modulemanager::ModulePtr<pulsar::modulebase::OneElectronIntegral> & mod,
                      const pulsar::system::BasisSet & bs);

std::vector<double>
FillTwoElectronVector(pulsar::modulemanager::ModulePtr<pulsar::modulebase::TwoElectronIntegral> & mod,
                      const pulsar::system::BasisSet & bs);

pulsar::math::IrrepSpinVectorD FindOccupations(size_t nelec);

pulsar::math::IrrepSpinMatrixD
FormDensity(const pulsar::math::IrrepSpinMatrixD & Cmat, const pulsar::math::IrrepSpinVectorD & occ);

double CalculateRMSDens(const pulsar::math::IrrepSpinMatrixD & m1,
                        const pulsar::math::IrrepSpinMatrixD & m2);

} // close namespace pulsarmethods

#endif

