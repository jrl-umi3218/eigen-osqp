/*
 * Copyright 2018-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Eigen_OSQP
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "eigen-osqp/OSQP.h"
#include "systems.h"

BOOST_FIXTURE_TEST_CASE(DENSE_MATRIX, QP1Dense)
{
  auto checkMatrix = [](const MatrixDense& matA, const MatrixDense& matConvA) {
    BOOST_REQUIRE_EQUAL(matA.rows(), matConvA.rows());
    BOOST_REQUIRE_EQUAL(matA.cols(), matConvA.cols());
    for (Eigen::Index i = 0; i < matA.rows(); ++i)
      for (Eigen::Index j = 0; j < matA.cols(); ++j)
        BOOST_REQUIRE_SMALL(std::abs(matA(i,j) - matConvA(i,j)), 1e-10);
  };

  Eigen::CSCMatrix cscA{};
  cscA.updateDefault(A);
  Eigen::MatrixXd convA = cscA.toDenseEigen();

  checkMatrix(A, convA);

  // Add block identity matrix underneath
  cscA.updateAndAddIdentity(A);
  Eigen::MatrixXd AId(A.rows() + nrvar, A.cols());
  AId.block(0, 0, A.rows(), A.cols()) = A;
  AId.block(A.rows(), 0, nrvar, nrvar).setIdentity();

  convA = cscA.toDenseEigen();

  checkMatrix(AId, convA);
}

BOOST_FIXTURE_TEST_CASE(SPARSE_MATRIX, QP1Sparse)
{
  auto checkMatrix = [](const MatrixSparse& matA, const MatrixSparse& matConvA) {
    BOOST_REQUIRE(matConvA.isCompressed());
    BOOST_REQUIRE_EQUAL(matA.rows(), matConvA.rows());
    BOOST_REQUIRE_EQUAL(matA.cols(), matConvA.cols());
    BOOST_REQUIRE_EQUAL(matA.outerSize(), matConvA.outerSize());
    for(int k = 0; k < matA.outerSize(); ++k)
    {
      for(MatrixSparse::InnerIterator itA(matA, k), itCA(matConvA, k); itA; ++itA, ++itCA)
      {
        BOOST_REQUIRE_EQUAL(itA.row(), itCA.row());
        BOOST_REQUIRE_EQUAL(itA.col(), itCA.col());
        BOOST_REQUIRE_SMALL(std::abs(itA.value() - itCA.value()), 1e-10);
      }
    }
  };

  Eigen::CSCMatrix cscA{};
  cscA.updateDefault(A);
  MatrixSparse convA = cscA.toSparseEigen();

  checkMatrix(A, convA);

  // Add block identity matrix underneath
  cscA.updateAndAddIdentity(A);

  Eigen::MatrixXd tmp(A.rows() + nrvar, A.cols());
  tmp.block(0, 0, A.rows(), A.cols()) = A;
  tmp.block(A.rows(), 0, nrvar, nrvar).setIdentity();
  MatrixSparse AId = tmp.sparseView();
  AId.makeCompressed();

  convA = cscA.toSparseEigen();

  checkMatrix(AId, convA);
}
