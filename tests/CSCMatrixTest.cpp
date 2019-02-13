#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Eigen_OSQP
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "eigen-osqp/OSQP.h"
#include "systems.h"

BOOST_FIXTURE_TEST_CASE(DENSE_MATRIX, QP1Dense)
{
  Eigen::CSCMatrix cscA(A);

  Eigen::MatrixXd convA = cscA.toDenseEigen();

  BOOST_REQUIRE_EQUAL(A.rows(), convA.rows());
  BOOST_REQUIRE_EQUAL(A.cols(), convA.cols());
  for (Eigen::Index i = 0; i < A.rows(); ++i)
    for (Eigen::Index j = 0; j < A.cols(); ++j)
      BOOST_REQUIRE_SMALL(std::abs(A(i,j) - convA(i,j)), 1e-10);

  // Add block identity matrix underneath
  cscA.updateAndAddIdentity(A);
  Eigen::MatrixXd AId(A.rows() + nrvar, A.cols());
  AId.block(0, 0, A.rows(), A.cols()) = A;
  AId.block(A.rows(), 0, nrvar, nrvar).setIdentity();

  convA = cscA.toDenseEigen();

  BOOST_REQUIRE_EQUAL(AId.rows(), convA.rows());
  BOOST_REQUIRE_EQUAL(AId.cols(), convA.cols());
  for (Eigen::Index i = 0; i < AId.rows(); ++i)
    for (Eigen::Index j = 0; j < AId.cols(); ++j)
      BOOST_REQUIRE_SMALL(std::abs(AId(i,j) - convA(i,j)), 1e-10);
}