#include "s21_matrix_oop.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <stdexcept>
#include <utility>
#include <utility>

S21Matrix::S21Matrix(int rows, int cols)
    : rows_ { rows }, cols_ { cols }, matrix_ {} {
  if (rows < 0) { throw std::range_error("Rows less than zero"); }
  if (cols < 0) { throw std::range_error("Cols less than zero"); }
  InitMatrix();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_ {other.rows_}, cols_ {other.cols_} {
  DeepCopy(other);
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_ { 0 }, cols_ { 0 }, matrix_ { nullptr } {
  std::swap(rows_, other.rows_);
  std::swap(cols_, other.cols_);
  std::swap(matrix_, other.matrix_);
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) {
    return *this;
  }
  DeepCopy(other);
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (&other == this) {
    return *this;
  }
  RemoveMatrix();
  rows_ = other.rows_;
  other.rows_ = 0;
  cols_ = other.cols_;
  other.cols_ = 0;
  matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  return *this;
}

S21Matrix::~S21Matrix() {
  RemoveMatrix();
}

double& S21Matrix::operator()(int row, int col) {
  if (row < 0 || row >= rows_) {
    throw std::range_error("Rows out of range");
  }
  if (col < 0 || col >= cols_) {
    throw std::range_error("Cols out of range");
  }
  return matrix_[row][col];
}

double S21Matrix::operator()(int row, int col) const {
  if (row < 0 || row >= rows_) {
    throw std::range_error("Rows out of range");
  }
  if (col < 0 || col >= cols_) {
    throw std::range_error("Cols out of range");
  }
  return matrix_[row][col];
}

void S21Matrix::set_rows(int rows) {
  ResizeMatrix(rows, cols_);
}

void S21Matrix::set_cols(int cols) {
  ResizeMatrix(rows_, cols);
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (int i { 0 }; i < rows_; ++i) {
    for (int j { 0 }; j < cols_; ++j) {
      if (std::fabs(matrix_[i][j] - other.matrix_[i][j]) > S21Matrix::kEps) {
          return false;
      }
    }
  }
  return true;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return this->EqMatrix(other);
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("SumMatrix matrix sizes are not equal");
  }
  for (int i { 0 }; i < rows_; ++i) {
    for (int j { 0 }; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::operator+=(const S21Matrix& other) {
  return this->SumMatrix(other);
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("SubMatrix matrix sizes are not equal");
  }
  for (int i { 0 }; i < rows_; ++i) {
    for (int j { 0 }; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::operator-=(const S21Matrix& other) {
  return this->SubMatrix(other);
}

void S21Matrix::MulNumber(const double num) {
  for (int i { 0 }; i < rows_; ++i) {
    for (int j { 0 }; j < cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::operator*=(const double num) {
  return this->MulNumber(num);
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::logic_error("In MulMatrix cols must be equal other rows");
  }
  S21Matrix tmpMatrix { rows_, other.cols_ };
  for (int i { 0 }; i < tmpMatrix.rows_; ++i) {
    for (int j { 0 }; j < tmpMatrix.cols_; ++j) {
      for (int k { 0 }; k < cols_; ++k) {
        tmpMatrix.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = std::move(tmpMatrix);
}

void S21Matrix::operator*=(const S21Matrix& other) {
  return this->MulMatrix(other);
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix tmpMatrix { cols_, rows_ };
  for (int i { 0 }; i < tmpMatrix.rows_; ++i) {
    for (int j { 0 }; j < tmpMatrix.cols_; ++j) {
        tmpMatrix.matrix_[i][j] = matrix_[j][i];
    }
  }
  return tmpMatrix;
}

S21Matrix S21Matrix::CalcCompements() const {
  if (rows_ != cols_) {
    throw std::logic_error("CalcCompements matrix should be square");
  }
  S21Matrix tmpMatrix { rows_, cols_ };
  if (rows_ == 1) {
    tmpMatrix(0, 0) = this->Determinant();
    return tmpMatrix;
  }
  for (int i { 0 }; i < rows_; ++i) {
    for (int j { 0 }; j < cols_; ++j) {
      tmpMatrix(i, j) = pow(-1, i + j) * (*this).Minor(i, j).Determinant();
    }
  }
  return tmpMatrix;
}

S21Matrix S21Matrix::Minor(int minor_i, int minor_j) const {
  if (minor_i < 0 || minor_j < 0) {
    throw std::range_error("Minor() arguments should be positive");
  }
  if (rows_ != cols_) {
    throw std::logic_error("Minor() matrix should be square");
  }
  S21Matrix tmpMatrix { *this };
  if (rows_ > 1) {
    tmpMatrix = S21Matrix { rows_ - 1, cols_ - 1 };
    for (int tmp_i { 0 }, i { 0 }; tmp_i < tmpMatrix.rows_; ++tmp_i, ++i) {
      if (tmp_i == minor_i) { ++i; }
      for (int tmp_j { 0 }, j { 0 }; tmp_j < tmpMatrix.cols_; ++tmp_j, ++j) {
        if (tmp_j == minor_j) { ++j; }
        tmpMatrix(tmp_i, tmp_j) = (*this)(i, j);
      }
    }
  }
  return tmpMatrix;
}

double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::logic_error("Can't calculate determinant for non-square matrix");
  }
  S21Matrix tmpMatrix{ *this };
  double determinant{ tmpMatrix.MakeMatrixTriangular() };
  for (int i = 0; i < rows_; ++i) {
    determinant *= tmpMatrix(i, i);
  }
  return determinant;
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (rows_ != cols_) {
    throw std::logic_error("InverseMatrix() matrix should be square");
  }
  if (std::fabs((*this).Determinant()) < S21Matrix::kEps) {
    throw std::logic_error("InverseMatrix() determinant == 0");
  }
  S21Matrix tmpMatrix { *this };
  if (rows_ == 1) {
    tmpMatrix(0, 0) = 1 / tmpMatrix(0, 0);
  } else {
    tmpMatrix = tmpMatrix.CalcCompements().Transpose();
    tmpMatrix.MulNumber(1 / (*this).Determinant());
  }
  return tmpMatrix;
}

void S21Matrix::InitMatrix() {
  matrix_ = new double*[rows_]{};
  for (int i { 0 }; i < rows_; ++i) {
    matrix_[i] = new double[cols_]{};
  }
}

void S21Matrix::RemoveMatrix() {
  if (matrix_ != nullptr) {  // To prevent nullptr dereference
    for (int i{ 0 }; i < rows_; ++i) {
      delete[] matrix_[i];
      matrix_[i] = nullptr;
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
}

void S21Matrix::CopyValuesFromMatrix(const S21Matrix& other) {
  for (int i { 0 }; i < (rows_ > other.rows_ ? other.rows_ : rows_); ++i) {
    for (int j { 0 }; j < (cols_ > other.cols_ ? other.cols_ : cols_); ++j) {
        matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

void S21Matrix::DeepCopy(const S21Matrix& other) {
  RemoveMatrix();
  rows_ = other.rows_;
  cols_ = other.cols_;
  InitMatrix();
  CopyValuesFromMatrix(other);
}

void S21Matrix::ResizeMatrix(int rows, int cols) {
  if (rows < 0) { throw std::range_error("Rows less than zero in resizing"); }
  if (cols < 0) { throw std::range_error("Cols less than zero in resizing"); }
  S21Matrix tmpMatrix { std::move(*this) };
  *this = S21Matrix{ rows, cols };
  CopyValuesFromMatrix(tmpMatrix);
}

double S21Matrix::MakeMatrixTriangular() {
  double determinant{ 1 };
  for (int i { 0 }; i < rows_ - 1; ++i) {
    int row{ i };
    while (row < rows_ - 1 &&
           std::fabs((*this)(row, i) < S21Matrix::kEps)) {
      ++row;
    }
    if (row == rows_ - 1 &&
        std::fabs((*this)(row, i) < S21Matrix::kEps)) {
      return 0;
    }
    if (row != i) {
      this->SwapRows(row, i);
      determinant = -determinant;
    }
    this->SubtractRows(i);
  }
  return determinant;
}

void S21Matrix::SwapRows(int row_1, int row_2) {
  std::swap(matrix_[row_1], matrix_[row_2]);
}

void S21Matrix::SubtractRows(int i) {
  for (int row{ i + 1 }; row < rows_; ++row) {
    double scale = (*this)(row, i) / (*this)(i, i);
    for (int j { i }; j < rows_; ++j) {
      (*this)(row, j) -= (*this)(i, j) * scale;
    }
  }
}

S21Matrix operator+(const S21Matrix& a, const S21Matrix& b) {
  S21Matrix tmpMatrix{ a };
  tmpMatrix += b;
  return tmpMatrix;
}

S21Matrix operator-(const S21Matrix& a, const S21Matrix& b) {
  return a + (b * -1);
}

S21Matrix operator*(const S21Matrix& a, double num) {
  S21Matrix tmpMatrix{ a };
  tmpMatrix *= num;
  return tmpMatrix;
}

S21Matrix operator*(double num, const S21Matrix& a) {
  return a * num;
}

S21Matrix operator*(const S21Matrix& a, const S21Matrix& b) {
  S21Matrix tmpMatrix{ a };
  tmpMatrix *= b;
  return tmpMatrix;
}

std::ostream& operator<<(std::ostream& out, const S21Matrix& matrix) {
  for (int i { 0 }; i < matrix.get_rows(); ++i) {
    for (int j { 0 }; j < matrix.get_cols(); ++j) {
        out << matrix(i, j) << ' ';
    }
    out << '\n';
  }
  return out;
}
