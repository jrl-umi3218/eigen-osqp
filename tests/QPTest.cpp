#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Eigen_OSQP
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <Eigen/Dense>

#include <eigen-osqp/OSQP.h>


struct QP1
{
  QP1()
  {
    nrvar = 6;
    nreq = 3;
    nrineq = 2;

    Q.resize(nrvar, nrvar);
    A.resize(nreq + nrineq, nrvar);

    C.resize(nrvar);
    AL.resize(nreq + nrineq);
    AU.resize(nreq + nrineq);
    XL.resize(nrvar);
    XU.resize(nrvar);
    X.resize(nrvar);

    auto inf = std::numeric_limits<double>::infinity();

    A << 1., -1., 1., 0., 3., 1.,
         -1., 0., -3., -4., 5., 6.,
         2., 5., 3., 0., 1., 0.,
         0., 1., 0., 1., 2., -1.,
         -1., 0., 2., 1., 1., 0.;
    AL << 1., 2., 3., -inf, -inf;
    AU << 1., 2., 3., -1., 2.5;

    //with  x between ci and cs:
    XL << -1000., -10000., 0., -1000., -1000.,-1000.;
    XU << 10000., 100., 1.5, 100., 100., 1000.;

    //and minimize 0.5*x'*Q*x + p'*x with
    C << 1., 2., 3., 4., 5., 6.;
    Q.setIdentity();

    X << 1.7975426, -0.3381487, 0.1633880, -4.9884023, 0.6054943, -3.1155623;
  }

  int nrvar, nreq, nrineq;
  Eigen::MatrixXd Q, A;
  Eigen::VectorXd C, AL, AU, XL, XU, X;
};


BOOST_AUTO_TEST_CASE(OSQP)
{
  QP1 qp1;

  Eigen::OSQP osqp{};

  osqp.problem(qp1.nrvar, qp1.nreq + qp1.nrineq);
  osqp.solve(qp1.Q, qp1.C,
    qp1.A, qp1.AL, qp1.AU,
    qp1.XL, qp1.XU);

  Eigen::VectorXd result = osqp.result();
  BOOST_CHECK_SMALL((osqp.result() - qp1.X).norm(), 1e-6);
}
