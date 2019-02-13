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

  osqp.problem(nrvar, nreq + nrineq);
  BOOST_REQUIRE(osqp.solve(Q, C, A, AL, AU, XL, XU));

  Eigen::VectorXd result = osqp.result();
  for (Eigen::Index i = 0; i < result.size(); ++i)
    BOOST_REQUIRE_SMALL(std::abs(result(i) - X(i)), 1e-6);
}
