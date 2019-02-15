#pragma once

#include "exportdecl.h"
#include "typedefs.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <osqp.h>

namespace Eigen
{

/*! \brief Helper class to convert Eigen::Matrix to csc struct expected by OSQP. */
struct EIGEN_OSQP_DLLAPI CSCMatrix
{
  /*! \brief Default constructor. Allocate memory for an empty matrix */
  CSCMatrix();

  /*! \brief Dense matrix to CSC conversion.
   * \param mat Dense matrix to convert.
   * \param doAddIdentity Ask to generate an identity matrix block under the converted matrix. Default is false.
   */
  CSCMatrix(const MatrixConstRef & mat, bool doAddIdentity = false);
  /*! \brief Sparse matrix to CSC conversion.
   * \param mat Sparse matrix to convert.
   * \param doAddIdentity Ask to generate an identity matrix block under the converted matrix. Default is false.
   */
  CSCMatrix(const MatrixCompressSparseConstRef & mat, bool doAddIdentity = false);
  /*! \brief Update current csc matrix to a new one.
   * \param mat Dense matrix to convert.
   * \param doAddIdentity Ask to generate an identity matrix block under the converted matrix. Default is false.
   */
  void update(const MatrixConstRef & mat, bool doAddIdentity = false);
  /*! \brief Update current csc matrix to a new one.
   * \param mat Sparse matrix to convert.
   * \param doAddIdentity Ask to generate an identity matrix block under the converted matrix. Default is false.
   */
  void update(const MatrixCompressSparseConstRef & mat, bool doAddIdentity = false);
  /*! \brief Get the csc matrix pointer */
  csc * matrix() noexcept { return &matrix_; }

  /*! \brief Convert csc matrix to Eigen dense matrix. */
  MatrixDense toDenseEigen() const;
  /*! \brief Convert csc matrix to Eigen sparse matrix. */
  MatrixSparse toSparseEigen() const;

private:
  void initParameters(const MatrixConstRef & mat, bool doAddIdentity = false);
  void initParameters(const MatrixCompressSparseConstRef & mat, bool doAddIdentity = false);
  void updateDefault(const MatrixConstRef & mat);
  void updateAndAddIdentity(const MatrixConstRef & mat);
  void updateDefault(const MatrixCompressSparseConstRef & mat);
  void updateAndAddIdentity(const MatrixCompressSparseConstRef & mat);

private:
  csc matrix_; /*!< OSQP sparse matrix representation */
  std::vector<c_int> p_; /*!< Vector of column index of the csc */
  std::vector<c_int> i_; /*!< Vector of row index of the csc */
  std::vector<c_float> x_; /*!< Vector of values of the csc */
};

} // namespace Eigen
