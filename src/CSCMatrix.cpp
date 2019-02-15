#include "eigen-osqp/CSCMatrix.h"
#include <numeric>

namespace Eigen
{

CSCMatrix::CSCMatrix()
{
  memset(&matrix_, 0, sizeof(csc));
}

CSCMatrix::CSCMatrix(const MatrixConstRef & mat, bool doAddIdentity)
: CSCMatrix()
{
  doAddIdentity ? updateAndAddIdentity(mat) : update(mat);
}

CSCMatrix::CSCMatrix(const MatrixCompressSparseConstRef & mat, bool doAddIdentity)
: CSCMatrix()
{
  doAddIdentity ? updateAndAddIdentity(mat) : update(mat);
}

void CSCMatrix::update(const MatrixConstRef & mat)
{
  initParameters(mat);

  Index nrRows = mat.rows();
  const auto* data = mat.data();
  auto iItr = i_.begin();

  // Generate column index
  std::generate(p_.begin(), p_.end(), [nrRows, n = 0]() mutable { return nrRows * (n++); });

  // Generate row index
  std::vector<int> rowIndex(nrRows);
  std::iota(rowIndex.begin(), rowIndex.end(), 0);
  for (Index j = 0; j < mat.cols(); ++j)
  {
    iItr = std::copy(rowIndex.cbegin(), rowIndex.cend(), iItr);
  }

  // Copy matrix value
  std::copy(data, data + mat.size(), x_.begin());
}

void CSCMatrix::updateAndAddIdentity(const MatrixConstRef & mat)
{
  initParameters(mat, true);

  Index i = 0;
  for(Index col_i = 0; col_i < mat.cols(); col_i++)
  {
    for(Index row_i = 0; row_i < mat.rows(); row_i++)
    {
      i_[i] = row_i;
      x_[i] = mat(row_i, col_i);
      i++;
    }
    p_[col_i] = col_i * (mat.rows() + 1);
    i_[i] = mat.rows() + col_i;
    x_[i] = 1.0;
    i++;
  }
}

void CSCMatrix::update(const MatrixCompressSparseConstRef & mat)
{
  initParameters(mat);

#if __cplusplus >= CPP17
  std::copy(std::execution::par_unseq, mat.outerIndexPtr(), mat.outerIndexPtr() + mat.outerSize(), p_.begin()); // Copy column start index of non-zeros
  std::copy(std::execution::par_unseq, mat.innerIndexPtr(), mat.innerIndexPtr() + mat.nonZeros(), i_.begin()); // Copy row index of non-zeros
  std::copy(std::execution::par_unseq, mat.valuePtr(), mat.valuePtr() + mat.nonZeros(), x_.begin()); // Copy matrix data
#else
  std::copy(mat.outerIndexPtr(), mat.outerIndexPtr() + mat.outerSize(), p_.begin()); // Copy column start index of non-zeros
  std::copy(mat.innerIndexPtr(), mat.innerIndexPtr() + mat.nonZeros(), i_.begin()); // Copy row index of non-zeros
  std::copy(mat.valuePtr(), mat.valuePtr() + mat.nonZeros(), x_.begin()); // Copy matrix data
#endif
}

void CSCMatrix::updateAndAddIdentity(const MatrixCompressSparseConstRef & mat)
{
  initParameters(mat, true);

  const auto* innerIndexIndexPtr = mat.innerIndexPtr();
  const auto* outerIndexIndexPtr = mat.outerIndexPtr();
  const auto* valueptr = mat.valuePtr();
  auto iItr = i_.begin();
  auto xItr = x_.begin();
  auto pItr = p_.begin();
  for (int k = 0; k < mat.outerSize(); ++k)
  {
    // Copy columns
    *(pItr++) = *(outerIndexIndexPtr++) + k; // A value from the identity matrix is added for each column

    int innerNNZs = *outerIndexIndexPtr - *(outerIndexIndexPtr - 1); // Get number of non-zeros

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

void CSCMatrix::initParameters(const MatrixConstRef & mat, bool doAddIdentity)
{
  Index nrVar = 0;
  if (doAddIdentity)
  {
    nrVar = mat.cols();
  }
  if(matrix_.nzmax != nrVar + mat.size())
  {
    matrix_.nzmax = nrVar + mat.size();
    matrix_.m = nrVar + mat.rows();
    matrix_.n = mat.cols();
    p_.resize(mat.cols() + 1);
    matrix_.p = p_.data();
    i_.resize(nrVar + mat.size());
    p_.back() = i_.size();
    matrix_.i = i_.data();
    x_.resize(nrVar + mat.size());
    matrix_.x = x_.data();
    matrix_.nz = -1;
  }
}

void CSCMatrix::initParameters(const MatrixCompressSparseConstRef & mat, bool doAddIdentity)
{
  Index nrVar = 0;
  if (doAddIdentity)
  {
    nrVar = mat.cols();
  }
  if(matrix_.nzmax != nrVar + mat.nonZeros())
  {
    matrix_.nzmax = nrVar + mat.nonZeros();
    matrix_.m = nrVar + mat.rows();
    matrix_.n = mat.cols();
    p_.resize(mat.cols() + 1);
    matrix_.p = p_.data();
    i_.resize(nrVar + mat.nonZeros());
    p_.back() = i_.size();
    matrix_.i = i_.data();
    x_.resize(nrVar + mat.nonZeros());
    matrix_.x = x_.data();
    matrix_.nz = -1;
  }
}

}
