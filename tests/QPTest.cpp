#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Eigen_OSQP
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <Eigen/Dense>

#include "eigen-osqp/OSQP.h"
#include "systems.h"

BOOST_FIXTURE_TEST_CASE(OSQP_DENSE, QP1Dense)
{
  Eigen::OSQP osqp{};

  // Test 1st solve function
  osqp.problem(nrvar, nreq + nrineq);
  BOOST_REQUIRE(osqp.solve(Q, C, A, AL, AU, XL, XU));

  Eigen::VectorXd result = osqp.result();
  for (Eigen::Index i = 0; i < result.size(); ++i)
    BOOST_REQUIRE_SMALL(std::abs(result(i) - X(i)), 1e-6);

  // Test 2nd solve function
  Eigen::MatrixXd A2(A.rows() + A.cols(), A.cols());
  A2.block(0, 0, A.rows(), A.cols()) = A;
  A2.block(A.rows(), 0, A.cols(), A2.cols()).setIdentity();
  Eigen::VectorXd AL2(AL.rows() + XL.rows());
  AL2.head(AL.rows()) = AL;
  AL2.tail(XL.rows()) = XL;
  Eigen::VectorXd AU2(AU.rows() + XU.rows());
  AU2.head(AL.rows()) = AU;
  AU2.tail(XL.rows()) = XU;

  BOOST_REQUIRE(osqp.solve(Q, C, A2, AL2, AU2));
  
  result = osqp.result();
  for (Eigen::Index i = 0; i < result.size(); ++i)
    BOOST_REQUIRE_SMALL(std::abs(result(i) - X(i)), 1e-6);
}

BOOST_FIXTURE_TEST_CASE(OSQP_SPARSE, QP1Sparse)
{
  Eigen::OSQP osqp{};

  // Test 1st solve function
  osqp.problem(nrvar, nreq + nrineq);
  BOOST_REQUIRE(osqp.solve(Q, C, A, AL, AU, XL, XU));

  Eigen::VectorXd result = osqp.result();
  for (Eigen::Index i = 0; i < result.size(); ++i)
    BOOST_REQUIRE_SMALL(std::abs(result(i) - X(i)), 1e-6);

  // Test 2nd solve function
  Eigen::MatrixXd tmp(A.rows() + A.cols(), A.cols());
  tmp.block(0, 0, A.rows(), A.cols()) = A;
  tmp.block(A.rows(), 0, A.cols(), tmp.cols()).setIdentity();
  Eigen::SparseMatrix<double> A2 = tmp.sparseView();
  Eigen::VectorXd AL2(AL.rows() + XL.rows());
  AL2.head(AL.rows()) = AL;
  AL2.tail(XL.rows()) = XL;
  Eigen::VectorXd AU2(AU.rows() + XU.rows());
  AU2.head(AL.rows()) = AU;
  AU2.tail(XL.rows()) = XU;

  BOOST_REQUIRE(osqp.solve(Q, C, A2, AL2, AU2));
  
  result = osqp.result();
  for (Eigen::Index i = 0; i < result.size(); ++i)
    BOOST_REQUIRE_SMALL(std::abs(result(i) - X(i)), 1e-6);
}