#include "eigen-osqp/CSCMatrix.h"
#include <numeric>

namespace Eigen
{

template <typename T>
Index triangular_size(T mat)
{
  assert(mat.rows() == mat.cols());

  Index size = 0;
  for (Index i = 1; i <= mat.rows(); ++i)
    size += i;

  return size;
}

CSCMatrix::CSCMatrix()
  : matrix_{0, 0, 0, nullptr, nullptr, nullptr, 0}
{
}

void CSCMatrix::updateTriangularDefault(const MatrixConstRef & mat)
{
  // Size of triangular matrix is sum_{i=1}^{N} i
  initParameters(mat.rows(), mat.cols(), triangular_size(mat), 0);

  Index i = 0;
  Eigen::Index cols = mat.cols();
  Eigen::Index rows = mat.rows();
  for(Index col = 0; col < cols; ++col) // Faster than using std::copy
  {
    p_[col] = i;
    for(Index row = 0; row < rows && row <= col; ++row)
    {
      i_[i] = row;
      x_[i] = mat(row, col);
      ++i;
    }
  }
}

void CSCMatrix::updateDefault(const MatrixConstRef & mat)
{
  initParameters(mat.rows(), mat.cols(), mat.size(), 0);

  Index i = 0;
  const auto* data = mat.data();
  Eigen::Index cols = mat.cols();
  Eigen::Index rows = mat.rows();
  for(Index col = 0; col < cols; ++col) // Faster than using std::copy
  {
    p_[col] = i;
    for(Index row = 0; row < rows; ++row)
    {
      i_[i] = row;
      x_[i] = *(data++);
      ++i;
    }
  }
}

void CSCMatrix::updateAndAddIdentity(const MatrixConstRef & mat)
{
  initParameters(mat.rows(), mat.cols(), mat.size(), mat.cols());

  Index i = 0;
  const auto* data = mat.data();
  Eigen::Index cols = mat.cols();
  Eigen::Index rows = mat.rows();
  for(Index col = 0; col < cols; ++col) // Faster than using std::copy
  {
    p_[col] = i;

    for(Index row = 0; row < rows; ++row)
    {
      i_[i] = row;
      x_[i] = *(data++);
      ++i;
    }

    i_[i] = rows + col;
    x_[i] = 1.0;
    ++i;
  }
}

void CSCMatrix::updateTriangularDefault(const MatrixCompressSparseConstRef & mat)
{
  // We can not initialize first since we don't know how many values are non-zero in the upper triangular part.
  // But we do know the maximum number of elements.
  Index maxValues = triangular_size(mat);
  i_.clear();
  x_.clear();
  p_.clear();
  i_.reserve(maxValues);
  x_.reserve(maxValues);
  p_.reserve(mat.cols() + 1);

  const auto* innerIndexIndexPtr = mat.innerIndexPtr();
  const auto* outerIndexIndexPtr = mat.outerIndexPtr();
  const auto* valueptr = mat.valuePtr();
  int innerNNZs;
  for (int k = 0; k < mat.outerSize(); ++k)
  {
    // Copy columns
    p_.push_back(*(outerIndexIndexPtr++));

    innerNNZs = *outerIndexIndexPtr - *(outerIndexIndexPtr - 1); // Get number of non-zeros

    // Copy rows
    for (int i = 0; i < innerNNZs; ++i, ++innerIndexIndexPtr, ++valueptr)
    {
      if (*innerIndexIndexPtr <= k) // Get only the upper part
      {
        i_.push_back(*innerIndexIndexPtr);
        x_.push_back(*valueptr);
      }
    }
  }

  // Initialize all parameters correctly. Inserted values in matrices won't be erased.
  initParameters(mat.rows(), mat.cols(), x_.size(), 0);
}

void CSCMatrix::updateDefault(const MatrixCompressSparseConstRef & mat)
{
  initParameters(mat.rows(), mat.cols(), mat.nonZeros(), 0);

  std::copy_n(mat.outerIndexPtr(), mat.outerSize(), p_.begin()); // Copy column start index of non-zeros
  std::copy_n(mat.innerIndexPtr(), mat.nonZeros(), i_.begin()); // Copy row index of non-zeros
  std::copy_n(mat.valuePtr(), mat.nonZeros(), x_.begin()); // Copy matrix data
}

void CSCMatrix::updateAndAddIdentity(const MatrixCompressSparseConstRef & mat)
{
  initParameters(mat.rows(), mat.cols(), mat.nonZeros(), mat.cols());

  const auto* innerIndexIndexPtr = mat.innerIndexPtr();
  const auto* outerIndexIndexPtr = mat.outerIndexPtr();
  const auto* valueptr = mat.valuePtr();
  auto iItr = i_.begin();
  auto xItr = x_.begin();
  auto pItr = p_.begin();
  int innerNNZs;
  for (int k = 0; k < mat.outerSize(); ++k)
  {
    // Copy columns
    *(pItr++) = *(outerIndexIndexPtr++) + k; // A value from the identity matrix is added for each column

    innerNNZs = *outerIndexIndexPtr - *(outerIndexIndexPtr - 1); // Get number of non-zeros

    // Copy rows
    iItr = std::copy_n(innerIndexIndexPtr, innerNNZs, iItr); // Copy rows
    innerIndexIndexPtr += innerNNZs;
    *(iItr++) = mat.rows() + k; // Add the identity matrix row

    // values
    xItr = std::copy_n(valueptr, innerNNZs, xItr); // Copy values
    valueptr += innerNNZs;
    *(xItr++) = 1.; // Add the identity matrix value
  }
}

MatrixDense CSCMatrix::toDenseEigen() const
{
  MatrixDense ret = MatrixDense::Zero(matrix_.m, matrix_.n);
  for(size_t i = 0; i < p_.size() - 1; ++i)
  {
    auto pi = p_[i];
    auto pi_next = p_[i+1];
    for(auto j = pi; j < pi_next; ++j)
    {
      ret(i_[j], i) = x_[j];
    }
  }
  return ret;
}

MatrixSparse CSCMatrix::toSparseEigen() const
{
  MatrixSparse ret(matrix_.m, matrix_.n);
  ret.reserve(matrix_.nzmax);
  for (Index j = 0; j < static_cast<Index>(p_.size()) - 1; ++j)
  {
    ret.startVec(j);
    for (Index i = p_[j]; i < p_[j + 1]; ++i)
    {
      ret.insertBack(i_[i], j) = x_[i];
    }
  }

  ret.finalize();
  return ret;
}

/*
 * Private methods
 */

void CSCMatrix::initParameters(Index rows, Index cols, Index newSize, Index nrIdVar)
{
  if(matrix_.nzmax != nrIdVar + newSize)
  {
    matrix_.nzmax = nrIdVar + newSize;
    matrix_.m = nrIdVar + rows;
    matrix_.n = cols;
    p_.resize(cols + 1);
    matrix_.p = p_.data();
    i_.resize(nrIdVar + newSize);
    p_.back() = i_.size();
    matrix_.i = i_.data();
    x_.resize(nrIdVar + newSize);
    matrix_.x = x_.data();
    matrix_.nz = -1;
  }
}

}
