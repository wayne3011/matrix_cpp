//
// Created by ackbarca on 08.02.24.
//

#ifndef CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
#define CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H

#define CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
#include <algorithm>
#include <cstring>
#include <iostream>
class S21Matrix {
 private:
  // Attributes
  int rows_ = 0, cols_ = 0;  // Rows and columns
  double** matrix_ =
      nullptr;  // Pointer to the memory where the matrix is allocated

 public:
  S21Matrix();  // Default constructor
  S21Matrix(int rows, int cols);

  S21Matrix(const S21Matrix& other);

  S21Matrix(S21Matrix&& other) noexcept;

  ~S21Matrix();  // Destructor

  // mutators
  void SetCols(int val);
  void SetRows(int val);
  // accessors
  [[nodiscard]] int GetCols() const { return cols_; }
  [[nodiscard]] int GetRows() const { return rows_; }

  // Functions

  void SumMatrix(const S21Matrix& other);

  bool EqMatrix(const S21Matrix& other);

  void SubMatrix(const S21Matrix& other);

  void MulNumber(double num);

  void MulMatrix(const S21Matrix& other);

  S21Matrix Transpose();

  S21Matrix CalcComplements();

  double Determinant();

  S21Matrix InverseMatrix();

  // operators
  S21Matrix& operator=(const S21Matrix& b);
  S21Matrix& operator+=(const S21Matrix& b);
  S21Matrix& operator-=(const S21Matrix& b);
  S21Matrix& operator*=(const S21Matrix& b);
  S21Matrix& operator*=(const double& b);

  double& operator()(int i, int j);
  double operator()(int i, int j) const;
  bool operator==(const S21Matrix& other) { return EqMatrix(other); }

 private:
  double CalcComplement(int i, int j);
  void ClearMatrix();
  void InitMatrix();
};
std::ostream& operator<<(std::ostream& os, const S21Matrix& d);

S21Matrix operator+(const S21Matrix& m1, const S21Matrix& m2);
S21Matrix operator-(const S21Matrix& m1, const S21Matrix& m2);
S21Matrix operator*(const S21Matrix& m1, const S21Matrix& m2);
S21Matrix operator*(const S21Matrix& m1, const double& number);
S21Matrix operator*(const double& number, const S21Matrix& m1);

#endif  // CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
