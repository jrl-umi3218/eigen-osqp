#include <type_traits>
#if __cplusplus >= cpp17
#include <execution>
#endif

namespace Eigen
{

template <typename TQ, typename TA>
bool OSQP::solve(const TQ & Q, const VectorConstRef & c, const TA & A, const VectorConstRef & AL,
                 const VectorConstRef & AU, const VectorConstRef & XL, const VectorConstRef & XU)
{
  // Check data
  assert(Q.rows() == Q.cols());
  assert(Q.rows() == c.rows());
  assert(A.rows() == AL.rows());
  assert(A.cols() == Q.rows());
  assert(AL.rows() == AU.rows());
  assert(XL.rows() == XU.rows());

  P_.update(Q);
  data_.P = P_.matrix();
  A_.updateAndAddIdentity(A);
  data_.A = A_.matrix();

  q_ = c;

#if __cplusplus >= cpp17
  // Copy lower bound
  std::copy_n(std::execution::par, AL.data(), AL.size(), bl_.data());
  std::copy_n(std::execution::par, XL.data(), XL.size(), bl_.data() + AL.size());
  // Copy upper bound
  std::copy_n(std::execution::par, AU.data(), AU.size(), bu_.data());
  std::copy_n(std::execution::par, XU.data(), XU.size(), bu_.data() + AU.size());
#else
  // Copy lower bound
  std::copy_n(AL.data(), AL.size(), bl_.data());
  std::copy_n(XL.data(), XL.size(), bl_.data() + AL.size());
  // Copy upper bound
  std::copy_n(AU.data(), AU.size(), bu_.data());
  std::copy_n(XU.data(), XU.size(), bu_.data() + AU.size());
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

if (ret)
{
#if __cplusplus >= cpp17
  std::copy_n(std::execution::par, workspace_->solution->x, result_.size(), result_.data());
#else
  std::copy_n(workspace_->solution->x, result_.size(), result_.data());
#endif
}

  return ret;
}

} // Eigen
