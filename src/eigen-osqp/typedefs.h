#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <osqp/glob_opts.h> // For c_float

// Define c++17 value of __cplusplus macro (for readability)
#define cpp17 201703L

// Must be Column-major !! (For now)
using MatrixConstRef = Eigen::Ref<const Eigen::MatrixXd>;
using VectorConstRef = Eigen::Ref<const Eigen::VectorXd>;
using MatrixSparse = Eigen::SparseMatrix<c_float>;
using MatrixCompressSparseConstRef = Eigen::Ref<const MatrixSparse, Eigen::StandardCompressedFormat>;
