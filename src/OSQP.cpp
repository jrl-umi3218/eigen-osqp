#include "eigen-osqp/OSQP.h"

void OSQPWorkspaceDeleter::operator()(OSQPWorkspace * ws) const
{
  osqp_cleanup(ws);
}

namespace Eigen
{

OSQP::OSQP()
: workspace_(nullptr)
{
  //FIXME Check later for better settings
  osqp_set_default_settings(&settings_);
  settings_.polish = 1;
  settings_.verbose = 0;
  memset(&data_, 0, sizeof(OSQPData));
}

OSQP::~OSQP()
{
}

bool OSQP::solve(const MatrixConstRef & Q,
                 const VectorConstRef & C,
                 const MatrixConstRef & A,
                 const VectorConstRef & AL,
                 const VectorConstRef & AU,
                 const VectorConstRef & XL,
                 const VectorConstRef & XU)
{
  assert(Q.rows() == Q.cols());
  assert(Q.rows() == C.rows());
  assert(A.rows() == AL.rows());
  assert(A.cols() == Q.rows());
  assert(AL.rows() == AU.rows());
  assert(XL.rows() == XU.rows());
  bool init_workspace = workspace_ == nullptr;
  // Update problem data
  if(data_.n != Q.rows())
  {
    data_.n = Q.rows();
    init_workspace = true;
  }
  if(data_.m != A.rows() + XL.rows())
  {
    data_.m = A.rows() + XL.rows();
    init_workspace = true;
  }
  P_.update(Q);
  data_.P = P_.matrix();
  A_.update(XL.rows(), A);
  data_.A = A_.matrix();
  if(init_workspace)
  {
    q_.resize(data_.n);
    bl_.resize(data_.m);
    bu_.resize(data_.m);
    result_.resize(data_.n);
  }
  memcpy(q_.data(), C.data(), C.size()*sizeof(double));
  data_.q = q_.data();
  memcpy(bl_.data(), XL.data(), XL.rows()*sizeof(double));
  memcpy(bu_.data(), XU.data(), XU.rows()*sizeof(double));
  memcpy(&(bl_.data()[XL.rows()]), AL.data(), AL.rows()*sizeof(double));
  memcpy(&(bu_.data()[XU.rows()]), AU.data(), AU.rows()*sizeof(double));
  data_.l = bl_.data();
  data_.u = bu_.data();
  // Initialize workspace_ if necessary
  if(init_workspace)
  {
    workspace_.reset(osqp_setup(&data_, &settings_));
  }
  // Solve
  bool ret = osqp_solve(workspace_.get()) >= 0;
  //if(ret)
  {
    memcpy(result_.data(), workspace_->solution->x, result_.size()*sizeof(double));
  }
  return ret;
}

}
