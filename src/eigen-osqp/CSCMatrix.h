#pragma once

#include <Eigen/Core>

#include <osqp/osqp.h>

using MatrixConstRef = Eigen::Ref<const Eigen::MatrixXd>;
using VectorConstRef = Eigen::Ref<const Eigen::VectorXd>;

namespace Eigen
{

/** Helper class to convert Eigen::Matrix to csc struct
 * expected by OSQP */
struct CSCMatrix
{
  /** Default (empty) matrix */
  CSCMatrix();

  /** Dense matrix to CSC conversion */
  CSCMatrix(const MatrixConstRef & mat);

  /** Block identity + dense matrix to CSC conversion */
  CSCMatrix(long long identity, const MatrixConstRef & mat);

  void update(const MatrixConstRef & mat);

  void update(long long identity, const MatrixConstRef & mat);

  csc * matrix() { return &matrix_; }

  /** For debugging */
  Eigen::MatrixXd toEigen() const;
private:
  csc matrix_;
  std::vector<c_int> p_;
  std::vector<c_int> i_;
  std::vector<c_float> x_;
  void update(const MatrixConstRef & mat, long long start);
};

} // namespace Eigen
