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
, data_{0, 0, nullptr, nullptr, nullptr, nullptr, nullptr}
{
  //FIXME Check later for better settings
  osqp_set_default_settings(&settings_);
  polish(true);
  verbose(false);
}

OSQP::~OSQP() noexcept
{
  if (data_.A)
  {
    delete data_.A;
  }
  if (data_.P)
  {
    delete data_.P;
  }
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
