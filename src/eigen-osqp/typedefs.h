#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <glob_opts.h> // For c_float

// Define c++17 value of __cplusplus macro (for readability)
#define cpp17 201703L

// Must be Column-major !! (For now)
using MatrixDense = Eigen::Matrix<c_float, Eigen::Dynamic, Eigen::Dynamic>;
using VectorDense = Eigen::Matrix<c_float, Eigen::Dynamic, 1>;
using MatrixConstRef = Eigen::Ref<const MatrixDense>;
using VectorConstRef = Eigen::Ref<const VectorDense>;
using MatrixSparse = Eigen::SparseMatrix<c_float>;
using MatrixCompressSparseConstRef = Eigen::Ref<const MatrixSparse, Eigen::StandardCompressedFormat>;
