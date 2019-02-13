#pragma once

#include "exportdecl.h"

#include <memory>

#include <Eigen/Core>

#include <eigen-osqp/CSCMatrix.h>

struct EIGEN_OSQP_DLLAPI OSQPWorkspaceDeleter
{
  void operator()(OSQPWorkspace * ptr) const;
};

namespace Eigen
{

class EIGEN_OSQP_DLLAPI OSQP
{
public:
  OSQP();

  ~OSQP();

  void problem(int nrVar, int nrConstr);

  const VectorDense & result() const { return result_; }

  template <typename TQ, typename TA>
  bool solve(const TQ & Q, const VectorConstRef & C, const TA & A, const VectorConstRef & AL, 
             const VectorConstRef & AU, const VectorConstRef & XL, const VectorConstRef & XU);

private:
  bool doInitWorkspace_;
  std::unique_ptr<OSQPWorkspace, OSQPWorkspaceDeleter> workspace_;
  /** Result */
  VectorDense result_;
  /** P matrix */
  CSCMatrix P_;
  /** A matrix */
  CSCMatrix A_;
  /** q: linear part of cost */
  VectorDense q_;
  /** Full lower bound (XL/AL) */
  VectorDense bl_;
  /** Full upper bound (XU/AU) */
  VectorDense bu_;
  /** Data of the solver */
  OSQPData data_;
  /** Settings */
  OSQPSettings settings_;
};

} // namespace Eigen

#include "OSQP.tpp"