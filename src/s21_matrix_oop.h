#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <exception>
#include <iostream>
#include <utility>

class S21Matrix {
 public:
  explicit S21Matrix(int rows = 3, int cols = 3);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other);
  ~S21Matrix();
  double& operator()(int row, int col);
  double operator()(int row, int col) const;
  inline int get_rows() const { return rows_; }
  inline int get_cols() const { return cols_; }
  void set_rows(int rows);
  void set_cols(int cols);
  bool EqMatrix(const S21Matrix& other) const;
  bool operator==(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void operator+=(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void MulNumber(const double num);
  void operator*=(const double num);
  void MulMatrix(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  S21Matrix Transpose() const;
  S21Matrix CalcCompements() const;
  S21Matrix Minor(int i, int j) const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;
  static constexpr double kEps { 1e-8 };

 private:
  void InitMatrix();
  void RemoveMatrix();
  void CopyValuesFromMatrix(const S21Matrix& other);
  void DeepCopy(const S21Matrix& other);
  void ResizeMatrix(int rows, int cols);
  double MakeMatrixTriangular();
  void SwapRows(int row_1, int row_2);
  void SubtractRows(int i);
  int rows_{}, cols_{};
  double** matrix_{};
};

S21Matrix operator+(const S21Matrix& a, const S21Matrix& b);
S21Matrix operator-(const S21Matrix& a, const S21Matrix& b);
S21Matrix operator*(const S21Matrix& a, double num);
S21Matrix operator*(double num, const S21Matrix& a);
S21Matrix operator*(const S21Matrix& a, const S21Matrix& b);

std::ostream& operator<<(std::ostream& out, const S21Matrix& matrix);
#endif  // SRC_S21_MATRIX_OOP_H_

