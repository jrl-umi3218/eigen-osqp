/*
 * Copyright 2018-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <type_traits>

namespace Eigen
{

template <typename TQ, typename TA>
bool OSQP::solve(const TQ & Q, const VectorConstRef & c, const TA & A, const VectorConstRef & AL,
                 const VectorConstRef & AU, const VectorConstRef & XL, const VectorConstRef & XU)
{
  static_assert(std::is_convertible<TQ, MatrixDense>::value || std::is_convertible<TQ, MatrixSparse>::value,
                "The type of Q should be at least convertible to an Eigen dense or sparse matrix of type c_float");
  static_assert(std::is_convertible<TA, MatrixDense>::value || std::is_convertible<TA, MatrixSparse>::value,
                "The type of A should be at least convertible to an Eigen dense or sparse matrix of type c_float");

  // Check data
  assert(Q.rows() == Q.cols());
  assert(Q.rows() == c.rows());
  assert(A.rows() == AL.rows());
  assert(A.cols() == Q.rows());
  assert(AL.rows() == AU.rows());
  assert(XL.rows() == XU.rows());

  P_.updateTriangularDefault(Q);
  data_.P = P_.matrix();
  A_.updateAndAddIdentity(A);
  data_.A = A_.matrix();

  auto blItr = bl_.begin();
  auto buItr = bu_.begin();
  // Copy linear part of cost
  std::copy_n(c.data(), c.size(), q_.begin());
  // Copy lower bound
  blItr = std::copy_n(AL.data(), AL.size(), blItr);
  std::copy_n(XL.data(), XL.size(), blItr);
  // Copy upper bound
  buItr = std::copy_n(AU.data(), AU.size(), buItr);
  std::copy_n(XU.data(), XU.size(), buItr);

  data_.q = q_.data();
  data_.l = bl_.data();
  data_.u = bu_.data();

  // Initialize workspace_ if necessary
  if(doInitWorkspace_)
  {
    OSQPWorkspace* tmpWork;
    doInitWorkspace_ = osqp_setup(&tmpWork, &data_, &settings_); // return 0 if success
    workspace_.reset(tmpWork);
    tmpWork = nullptr;
  }

  // Solve
  bool ret = osqp_solve(workspace_.get()) >= 0;

  if (ret)
  {
    std::copy_n(workspace_->solution->x, result_.size(), result_.data());
  }

  return ret;
}

template <typename TQ, typename TA>
bool OSQP::solve(const TQ & Q, const VectorConstRef & c, const TA & A, const VectorConstRef & AL, const VectorConstRef & AU)
{
  static_assert(std::is_convertible<TQ, MatrixDense>::value || std::is_convertible<TQ, MatrixSparse>::value,
                "The type of Q should be at least convertible to an Eigen dense or sparse matrix of type c_float");
  static_assert(std::is_convertible<TA, MatrixDense>::value || std::is_convertible<TA, MatrixSparse>::value,
                "The type of A should be at least convertible to an Eigen dense or sparse matrix of type c_float");

  // Check data
  assert(Q.rows() == Q.cols());
  assert(Q.rows() == c.rows());
  assert(A.rows() == AL.rows());
  assert(A.cols() == Q.rows());

  P_.updateTriangularDefault(Q);
  data_.P = P_.matrix();
  A_.updateDefault(A);
  data_.A = A_.matrix();

  // Copy linear part of cost
  std::copy_n(c.data(), c.size(), q_.begin());
  // Copy lower bound
  std::copy_n(AL.data(), AL.size(), bl_.begin());
  // Copy upper bound
  std::copy_n(AU.data(), AU.size(), bu_.begin());

  data_.q = q_.data();
  data_.l = bl_.data();
  data_.u = bu_.data();

  // Initialize workspace_ if necessary
  if(doInitWorkspace_)
  {
    OSQPWorkspace* tmpWork;
    doInitWorkspace_ = osqp_setup(&tmpWork, &data_, &settings_); // return 0 if success
    workspace_.reset(tmpWork);
    tmpWork = nullptr;
  }

  // Solve
  bool ret = osqp_solve(workspace_.get()) >= 0;

  if (ret)
  {
    std::copy_n(workspace_->solution->x, result_.size(), result_.data());
  }

  return ret;
}

template <typename TQ>
bool OSQP::solve(const TQ & Q, const VectorConstRef & c)
{
  static_assert(std::is_convertible<TQ, MatrixDense>::value || std::is_convertible<TQ, MatrixSparse>::value,
                "The type of Q should be at least convertible to an Eigen dense or sparse matrix of type c_float");

  // Check data
  assert(Q.rows() == Q.cols());
  assert(Q.rows() == c.rows());

  P_.updateTriangularDefault(Q);
  data_.P = P_.matrix();

  // Copy linear part of cost
  std::copy_n(c.data(), c.size(), q_.begin());

  data_.q = q_.data();

  // Initialize workspace_ if necessary
  if(doInitWorkspace_)
  {
    OSQPWorkspace* tmpWork;
    doInitWorkspace_ = osqp_setup(&tmpWork, &data_, &settings_); // return 0 if success
    workspace_.reset(tmpWork);
    tmpWork = nullptr;
  }

  // Solve
  bool ret = osqp_solve(workspace_.get()) >= 0;

  if (ret)
  {
    std::copy_n(workspace_->solution->x, result_.size(), result_.data());
  }

  return ret;
}

} // Eigen
