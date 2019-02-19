#include "eigen-osqp/OSQP.h"

void OSQPWorkspaceDeleter::operator()(OSQPWorkspace * ws) const
{
  osqp_cleanup(ws);
}

namespace Eigen
{

OSQP::OSQP()
: doInitWorkspace_(true)
, workspace_(nullptr)
{
  //FIXME Check later for better settings
  osqp_set_default_settings(&settings_);
  polish(true);
  verbose(false);
  memset(&data_, 0, sizeof(OSQPData));
}

void OSQP::problem(int nrVar, int nrConstr)
{
  if(data_.n != nrVar || data_.m != nrConstr)
  {
    data_.n = nrVar;
    data_.m = nrVar + nrConstr;
    q_.resize(data_.n);
    bl_.resize(data_.m);
    bu_.resize(data_.m);
    result_.resize(data_.n);
    doInitWorkspace_ = true;
  }
}

}
