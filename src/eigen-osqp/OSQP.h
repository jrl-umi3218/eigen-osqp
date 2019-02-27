#pragma once

#include "exportdecl.h"

#include <memory>

#include <Eigen/Core>

#include <eigen-osqp/CSCMatrix.h>

/*! \brief Deleter object for OSQPWorkspace. */
struct EIGEN_OSQP_DLLAPI OSQPWorkspaceDeleter
{
  void operator()(OSQPWorkspace * ptr) const;
};

namespace Eigen
{

/*! \brief Eigen wrapper class of osqp.
 * The class wraps the c-library osqp for Eigen matrix.
 * It accepts both dense and sparse matrix.
 * To compute a QP, just set the number of variables and solve.
 * 
 * \note Dense matrix are converted to an osqp sparse representation.
 * All zeros in it won't be deleted away.
 */
class EIGEN_OSQP_DLLAPI OSQP
{
public:
  /*! \brief Default constructor. */
  OSQP();
  /*! \brief Default destructor. */
  ~OSQP() noexcept = default;

  // Access to settings.
  /*! \brief ADMM step rho. Default is 0.1. */
  void admmStepRho(c_float rho) noexcept { settings_.rho = rho; }
  c_float admmStepRho() const noexcept { return settings_.rho; }

  /*! \brief ADMM step sigma. Default is 1e-6. */
  void admmStepSigma(c_float sigma) noexcept { settings_.sigma = sigma; }
  c_float admmStepSigma() const noexcept { return settings_.sigma; }

  /*! \brief Heuristic data scaling iterations. If 0, scaling disabled. Default is 10. */
  void scalingIter(c_int scaling) noexcept { settings_.scaling = scaling; }
  c_int scalingIter() const noexcept { return settings_.scaling; }

  /*! \brief Maximum number of iteration. Default is 4000. */
  void maxIter(c_int maxIter) noexcept { settings_.max_iter = maxIter; }
  c_int maxIter() const noexcept { return settings_.max_iter; }

  /*! \brief Absolute convergence tolerance. Default is 1e-3. */
  void absConvergenceTol(c_float tol) noexcept { settings_.eps_abs = tol; }
  c_float absConvergenceTol() const noexcept { return settings_.eps_abs; }

  /*! \brief Relative convergence tolerance. Default is 1e-3. */
  void relConvergenceTol(c_float tol) noexcept { settings_.eps_rel = tol; }
  c_float relConvergenceTol() const noexcept { return settings_.eps_rel; }

  /*! \brief Primal infeasibility tolerance. Default is 1e-4. */
  void primalInfeasibilityTol(c_float tol) noexcept { settings_.eps_prim_inf = tol; }
  c_float primalInfeasibilityTol() const noexcept { return settings_.eps_prim_inf; }

  /*! \brief Dual infeasibility tolerance. Default is 1e-4. */
  void dualInfeasibilityTol(c_float tol) noexcept { settings_.eps_dual_inf = tol; }
  c_float dualInfeasibilityTol() const noexcept { return settings_.eps_dual_inf; }

  /*! \brief Relaxation parameter. Default is 1.6. */
  void relaxationParam(c_float alpha) noexcept { settings_.alpha = alpha; }
  c_float relaxationParam() const noexcept { return settings_.alpha; }

  /*! \brief Linear solver. Default is QDLDL_SOLVER. */
  void linearSystemSolver(linsys_solver_type solver) noexcept { settings_.linsys_solver = solver; }
  linsys_solver_type linearSystemSolver() const noexcept { return settings_.linsys_solver; }

  /*! \brief Use scale termination criteria. Default is false. */
  void scaleTermination(bool doScaleTermination) noexcept { settings_.scaled_termination = (doScaleTermination ? 1 : 0); }
  bool scaleTermination() const noexcept { return static_cast<bool>(settings_.scaled_termination); }

  /*! \brief check termination interval. If 0, termination checking is disabled. Default is 25. */
  void checkTermination(c_int checkValue) noexcept { settings_.check_termination = checkValue; }
  c_int checkTermination() const noexcept { return settings_.check_termination; }

  /*! \brief Warm start. Default is true. */
  void warmStart(bool ws) noexcept { settings_.warm_start = (ws ? 1 : 0); }
  bool warmStart() const noexcept { return static_cast<bool>(settings_.warm_start); }

  /*! \brief Polish ADMM solution. */
  void polish(bool doPolish) noexcept { settings_.polish = (doPolish ? 1 : 0); }
  bool polish() const noexcept { return static_cast<bool>(settings_.polish); }

  /*! \brief Regularization parameter for polish. */
  void polishDelta(c_float delta) noexcept { settings_.delta = delta; }
  c_float polishDelta() const noexcept { return settings_.delta; }

  /*! \brief Iterative refinement steps in polish. */
  void polishRefineIter(c_int iter) noexcept { settings_.polish_refine_iter = iter; }
  c_int polishRefineIter() const noexcept { return settings_.polish_refine_iter; }

  /*! \brief Write osqp progress. */
  void verbose(bool doPrint) noexcept { settings_.verbose = (doPrint ? 1 : 0); }
  bool verbose() const noexcept { return static_cast<bool>(settings_.verbose); }

  /*! \brief Set up the problem parameters.
   * \param nrVar Number of variable
   * \param nrConstr Number of constraints
   */
  void problem(int nrVar, int nrConstr);
  /*! \brief Return the qp result.
   * \return Result (uninitialize vector if the qp fails)
   */
  const VectorDense & result() const noexcept { return result_; }
  /*! \brief Solve the given problem.
   * The function problem(int, int) needs to be called before calling this function.
   * \tparam TQ Matrix type. Either a dense or compressed sparse matrix (or a reference to it).
   * \tparam TA Matrix type. Either a dense or compressed sparse matrix (or a reference to it).
   * \param Q Quadratic part of the cost.
   * \param c Linear part of the cost.
   * \param A Constraint matrix.
   * \param AL Lower bound of the constraints.
   * \param AU Upper bound of the constraints.
   * \param XL Lower bound of the variables.
   * \param XU Upper bound of the variables.
   */
  template <typename TQ, typename TA>
  bool solve(const TQ & Q, const VectorConstRef & c, const TA & A, const VectorConstRef & AL, 
             const VectorConstRef & AU, const VectorConstRef & XL, const VectorConstRef & XU);
  /*! \brief Solve the given problem.
   * The function problem(int, int) needs to be called before calling this function.
   * \tparam TQ Matrix type. Either a dense or compressed sparse matrix (or a reference to it).
   * \tparam TA Matrix type. Either a dense or compressed sparse matrix (or a reference to it).
   * \param Q Quadratic part of the cost.
   * \param c Linear part of the cost.
   * \param A Constraint matrix.
   * \param AL Lower bound of the constraints.
   * \param AU Upper bound of the constraints.
   */
  template <typename TQ, typename TA>
  bool solve(const TQ & Q, const VectorConstRef & c, const TA & A, const VectorConstRef & AL, const VectorConstRef & AU);
  /*! \brief Solve unconstrained problem.
   * The function problem(int, int) needs to be called before calling this function.
   * \tparam TQ Matrix type. Either a dense or compressed sparse matrix (or a reference to it).
   * \param Q Quadratic part of the cost.
   * \param c Linear part of the cost.
   */
  template <typename TQ>
  bool solve(const TQ & Q, const VectorConstRef & c);

  // Results info
  /*! \brief Get number of iterations taken. */
  c_int iter() const noexcept { return workspace_->info->iter; }
  /*! \brief Print solver status. */
  void inform(std::ostream& os) const noexcept { os << workspace_->info->status << "\n"; }
  /*! \brief Get solver status. */
  c_int status() const noexcept { return workspace_->info->status_val; }
  /*! \brief Get cost result. */
  c_float costResult() const noexcept { return workspace_->info->obj_val; }
  /*! \brief Get norm of primal residual. */
  c_float primalResidualNorm() const noexcept { return workspace_->info->pri_res; }
  /*! \brief Get norm of dual residual. */
  c_float dualResidualNorm() const noexcept { return workspace_->info->dua_res; }

private:
  bool doInitWorkspace_; /*!< True if the workspace needs to be initialized */
  std::unique_ptr<OSQPWorkspace, OSQPWorkspaceDeleter> workspace_; /*!< Workspace of osqp */
  VectorDense result_; /*!< Result of optimization (Uninitialized if qp fails) */
  CSCMatrix P_; /*!< Quadratic part of the cost */
  CSCMatrix A_; /*!< Constraint matrix */
  std::vector<c_float> q_; /*!< Linear part of the cost */
  std::vector<c_float> bl_; /*!< Lower bound constraint */
  std::vector<c_float> bu_; /*!< Upper bound constraint */
  OSQPData data_; /*!< Data of the solver */
  OSQPSettings settings_; /*!< Solver settings */
};

} // namespace Eigen

#include "OSQP.hpp"