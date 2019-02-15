#include <type_traits>
#if __cplusplus >= CPP17
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

  auto blItr = bl_.begin();
  auto buItr = bu_.begin();
#if __cplusplus >= CPP17
  // Copy linear part of cost
  std::copy_n(std::execution::par_unseq, c.data(), c.size(), q_.begin());
  // Copy lower bound
  blItr = std::copy_n(std::execution::par_unseq, AL.data(), AL.size(), blItr);
  std::copy_n(std::execution::par_unseq, XL.data(), XL.size(), blItr);
  // Copy upper bound
  buItr = std::copy_n(std::execution::par_unseq, AU.data(), AU.size(), buItr);
  std::copy_n(std::execution::par_unseq, XU.data(), XU.size(), buItr);
#else
  // Copy lower bound
  blItr = std::copy_n(AL.data(), AL.size(), blItr);
  std::copy_n(XL.data(), XL.size(), blItr);
  // Copy upper bound
  buItr = std::copy_n(AU.data(), AU.size(), buItr);
  std::copy_n(XU.data(), XU.size(), buItr);
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
#if __cplusplus >= CPP17
  std::copy_n(std::execution::par_unseq, workspace_->solution->x, result_.size(), result_.data());
#else
  std::copy_n(workspace_->solution->x, result_.size(), result_.data());
#endif
}

  return ret;
}

} // Eigen
