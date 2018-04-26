#include "eigen-osqp/CSCMatrix.h"

namespace Eigen
{

CSCMatrix::CSCMatrix()
{
  memset(&matrix_, 0, sizeof(csc));
}

CSCMatrix::CSCMatrix(const MatrixConstRef & mat)
: CSCMatrix()
{
  update(mat);
}

CSCMatrix::CSCMatrix(long long identity, const MatrixConstRef & mat)
: CSCMatrix()
{
  update(identity, mat);
}

void CSCMatrix::update(const MatrixConstRef & mat)
{
  update(mat, 0);
}

void CSCMatrix::update(long long identity, const MatrixConstRef & mat)
{
  update(mat, identity);
}

void CSCMatrix::update(const MatrixConstRef & mat, long long identity)
{
  if(matrix_.nzmax != mat.size() + identity)
  {
    matrix_.nzmax = mat.size() + identity;
    matrix_.m = identity + mat.rows();
    matrix_.n = mat.cols();
    p_.resize(mat.cols() + 1);
    matrix_.p = p_.data();
    i_.resize(identity + mat.size());
    p_.back() = i_.size();
    matrix_.i = i_.data();
    x_.resize(identity + mat.size());
    matrix_.x = x_.data();
    matrix_.nz = -1;
  }
  size_t i = 0;
  for(long long col_i = 0; col_i < mat.cols(); col_i++)
  {
    if(identity > col_i)
    {
      p_[col_i] = col_i * (mat.rows() + 1);
      i_[i] = col_i;
      x_[i] = 1.0;
      i++;
    }
    else
    {
      p_[col_i] = col_i * mat.rows();
    }
    for(long long row_i = 0; row_i < mat.rows(); row_i++)
    {
      i_[i] = row_i + identity;
      x_[i] = mat(row_i, col_i);
      i++;
    }
  }
}

Eigen::MatrixXd CSCMatrix::toEigen() const
{
  MatrixXd ret = Eigen::MatrixXd::Zero(matrix_.m, matrix_.n);
  for(size_t i = 0; i < p_.size() - 1; ++i)
  {
    auto pi = p_[i];
    auto pi_next = p_[i+1];
    for(size_t j = pi; j < pi_next; ++j)
    {
      ret(i_[j], i) = x_[j];
    }
  }
  return ret;
}

}
