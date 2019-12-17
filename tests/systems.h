#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

struct QP1Dense
{
  QP1Dense()
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

struct QP1Sparse
{
  QP1Sparse()
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

    std::vector<Eigen::Triplet<double>> nz = {
        {0, 0, 1.},  {0, 1, -1.},  {0, 2, 1.},               {0, 4, 3.}, {0, 5, 1.},
        {1, 0, -1.},               {1, 2, -3.}, {1, 3, -4.}, {1, 4, 5.}, {1, 5, 6.},
        {2, 0, 2.},   {2, 1, 5.},  {2, 2, 3.},               {2, 4, 1.},
                      {3, 1, 1.},               {3, 3, 1.},  {3, 4, 2.}, {3, 5, -1.},
        {4, 0, -1.},               {4, 2, 2.},  {4, 3, 1.},  {4, 4, 1.},
    };

    A.setFromTriplets(nz.cbegin(), nz.cend());
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
  Eigen::SparseMatrix<double> Q, A;
  Eigen::VectorXd C, AL, AU, XL, XU, X;
};
