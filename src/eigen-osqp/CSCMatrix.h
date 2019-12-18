/*
 * Copyright 2018-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

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

  /*! \brief Update current csc matrix to a new one. Only the upper triangular part is considered.
   * \note It is up to the user to check that the matrix is symmetrical.
   * \param mat Dense matrix to convert.
   */
  void updateTriangularDefault(const MatrixConstRef & mat);
  /*! \brief Update current csc matrix to a new one.
   * \param mat Dense matrix to convert.
   */
  void updateDefault(const MatrixConstRef & mat);
  /*! \brief Update current csc matrix to a new one and add an identity matrix beneath it.
   * \param mat Dense matrix to convert.
   */
  void updateAndAddIdentity(const MatrixConstRef & mat);
  /*! \brief Update current csc matrix to a new one. Only the upper triangular part is considered.
   * \note It is up to the user to check that the matrix is symmetrical.
   * \param mat Dense matrix to convert.
   */
  void updateTriangularDefault(const MatrixCompressSparseConstRef & mat);
  /*! \brief Update current csc matrix to a new one.
   * \param mat Compress sparse matrix to convert.
   */
  void updateDefault(const MatrixCompressSparseConstRef & mat);
  /*! \brief Update current csc matrix to a new one and add an identity matrix beneath it.
   * \param mat Compress sparse matrix to convert.
   */
  void updateAndAddIdentity(const MatrixCompressSparseConstRef & mat);

  /*! \brief Get the csc matrix pointer */
  csc * matrix() noexcept { return &matrix_; }

  /*! \brief Convert csc matrix to Eigen dense matrix. */
  MatrixDense toDenseEigen() const;
  /*! \brief Convert csc matrix to Eigen sparse matrix. */
  MatrixSparse toSparseEigen() const;

private:
  void initParameters(Index rows, Index cols, Index newSize, Index nrIdVar);

private:
  csc matrix_; /*!< OSQP sparse matrix representation */
  std::vector<c_int> p_; /*!< Vector of column index of the csc */
  std::vector<c_int> i_; /*!< Vector of row index of the csc */
  std::vector<c_float> x_; /*!< Vector of values of the csc */
};

} // namespace Eigen
