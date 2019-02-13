#pragma once

#include "typedefs.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <osqp/osqp.h>

namespace Eigen
{

/** Helper class to convert Eigen::Matrix to csc struct
 * expected by OSQP */
struct CSCMatrix
{
  /** Default (empty) matrix */
  CSCMatrix();

  /** Dense matrix to CSC conversion */
  CSCMatrix(const MatrixConstRef & mat, bool doAddIdentity = false);

  /** Sparse matrix to CSC conversion */
  CSCMatrix(const MatrixCompressSparseConstRef & mat, bool doAddIdentity = false);

  void update(const MatrixConstRef & mat);
  void updateAndAddIdentity(const MatrixConstRef & mat);

  void update(const MatrixCompressSparseConstRef & mat);
  void updateAndAddIdentity(const MatrixCompressSparseConstRef & mat);

  csc * matrix() { return &matrix_; }

  /** For debugging */
  MatrixXd toDenseEigen() const;
  MatrixSparse toSparseEigen() const;

private:
  void initParameters(const MatrixConstRef & mat, bool doAddIdentity = false);
  void initParameters(const MatrixCompressSparseConstRef & mat, bool doAddIdentity = false);

private:
  csc matrix_;
  std::vector<c_int> p_;
  std::vector<c_int> i_;
  std::vector<c_float> x_;
};

} // namespace Eigen
