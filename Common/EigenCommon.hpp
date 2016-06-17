#ifndef PULSAR_GUARD_COMMON__EIGENCOMMON_HPP_
#define PULSAR_GUARD_COMMON__EIGENCOMMON_HPP_

#include <pulsar/math/BlockByIrrepSpin.hpp>
#include <pulsar/math/TensorImpl.hpp>

#include <Eigen/Dense>


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


std::shared_ptr<const Eigen::MatrixXd>
convert_to_eigen(const std::shared_ptr<const pulsar::math::TensorImpl<2, double>> & ten);

std::shared_ptr<const Eigen::VectorXd>
convert_to_eigen(const std::shared_ptr<const pulsar::math::TensorImpl<1, double>> & ten);


typedef pulsar::math::BlockByIrrepSpin<Eigen::MatrixXd> BlockedEigenMatrix;
typedef pulsar::math::BlockByIrrepSpin<Eigen::VectorXd> BlockedEigenVector;


#endif

