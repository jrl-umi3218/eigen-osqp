#include <type_traits>
#include <execution>

namespace Eigen
{

template <typename TQ, typename TA>
bool OSQP::solve(const TQ & Q, const VectorConstRef & c, const TA & A, const VectorConstRef & AL,
                 const VectorConstRef & AU, const VectorConstRef & XL, const VectorConstRef & XU)
{
  // Check data
  assert(Q.rows() == Q.cols());
  assert(Q.rows() == C.rows());
  assert(A.rows() == AL.rows());
  assert(A.cols() == Q.rows());
  assert(AL.rows() == AU.rows());
  assert(XL.rows() == XU.rows());

  P_.update(Q);
  data_.P = P_.matrix();
  A_.updateAndAddIdentity(A);
  data_.A = A_.matrix();

#if __cplusplus >= cpp17
  // Copy q
  std::copy(std::execution::par, c.data(), c.data() + c.size(), q_.begin());
  // Copy lower bound
  std::copy(std::execution::par, AL.data(), AL.data() + AL.size(), bl_.begin());
  std::copy(std::execution::par, XL.data(), XL.data() + XL.size(), bl_.begin() + AL.size());
  // Copy upper bound
  std::copy(std::execution::par, AU.data(), AU.data() + AU.size(), bu_.begin());
  std::copy(std::execution::par, XU.data(), XU.data() + XU.size(), bu_.begin() + AU.size());
#else
  // Copy q
  std::copy(c.data(), c.data() + c.size(), q_.begin());
  // Copy lower bound
  std::copy(AL.data(), AL.data() + AL.size(), bl_.begin());
  std::copy(XL.data(), XL.data() + XL.size(), bl_.begin() + AL.size());
  // Copy upper bound
  std::copy(AU.data(), AU.data() + AU.size(), bu_.begin());
  std::copy(XU.data(), XU.data() + XU.size(), bu_.begin() + AU.size());
#endif

  data_.q = q_.data();
  data_.l = bl_.data();
  data_.u = bu_.data();

  // Initialize workspace_ if necessary
  if(doInitWorkspace_)
  {
    workspace_.reset(osqp_setup(&data_, &settings_));
    doInitWorkspace_ = false;
  }

  // Solve
  bool ret = osqp_solve(workspace_.get()) >= 0;

#if __cplusplus >= cpp17
  std::copy(std::execution::par, workspace_->solution->x, workspace_->solution->x + result_.size(), result_.data());
#else
  std::copy(workspace_->solution->x, workspace_->solution->x + result_.size(), result_.data());
#endif

  return ret;
}

} // Eigen
