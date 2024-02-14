//
// Created by ackbarca on 08.02.24.
//

#include "s21_matrix_oop.h"

#include <cmath>
#include <stdexcept>
//      CONSTRUCTORS && DESTRUCTORS
S21Matrix::S21Matrix() : rows_(3), cols_(3) { InitMatrix(); }

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 1 || cols < 1) throw std::out_of_range("Incorrect matrix size!");
  InitMatrix();
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < other.rows_; ++i) {
    std::copy(other.matrix_[i], other.matrix_[i] + other.cols_, matrix_[i]);
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() { ClearMatrix(); }

//      ACCESSORS && MUTATORS

void S21Matrix::SetCols(int val) {
  if (cols_ == val) return;
  if (val < 1) throw std::out_of_range("Invalid columns count!");
  int cpy_length = val > cols_ ? cols_ : val;
  for (int i = 0; i < rows_; ++i) {
    double *last_row = matrix_[i];
    auto new_row = new double[val];
    std::copy(last_row, last_row + cpy_length, new_row);
    matrix_[i] = new_row;
    delete[] last_row;
  }
  cols_ = val;
}

void S21Matrix::SetRows(int val) {
  if (val == rows_) return;
  if (val < 1) throw std::out_of_range("Invalid rows count!");
  int cpy_length = val > rows_ ? rows_ : val;
  auto new_matrix = new double *[val];
  std::copy(matrix_, matrix_ + cpy_length, new_matrix);

  if (val > rows_) {
    for (int i = cpy_length; i < val; ++i) {
      new_matrix[i] = new double[cols_];
    }
  }
  if (val < rows_) {
    for (int i = val; i < rows_; ++i) {
      delete[] matrix_[i];
    }
  }
  delete[] matrix_;
  matrix_ = new_matrix;
  rows_ = val;
}

//      PUBLIC METHODS

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  if (other.rows_ != rows_ || other.cols_ != cols_) return false;

  for (int i = 0; i < other.rows_; ++i) {
    for (int j = 0; j < other.cols_; ++j) {
      if (std::fabs(matrix_[i][j] - other.matrix_[i][j]) >=
          std::numeric_limits<double>::epsilon())
        return false;
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (cols_ != other.cols_ || rows_ != other.rows_)
    throw std::invalid_argument(
        "The size of the other matrix does not match the current one!");

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) matrix_[i][j] += other.matrix_[i][j];
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (cols_ != other.cols_ || rows_ != other.rows_)
    throw std::invalid_argument(
        "The size of the other matrix does not match the current one!");

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) matrix_[i][j] -= other.matrix_[i][j];
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_)
    throw std::invalid_argument(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix.");

  S21Matrix thisCopy = S21Matrix(*this);

  SetCols(other.cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = 0;
      for (int k = 0; k < thisCopy.cols_; k++) {
        matrix_[i][j] += (thisCopy.matrix_[i][k] * other.matrix_[k][j]);
      }
    }
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result = S21Matrix(cols_, rows_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  double result = 1;
  if (rows_ != cols_) throw std::logic_error("The matrix is not a square!");
  int size = rows_;

  S21Matrix temp_matrix = S21Matrix(*this);

  for (int i = 0; i < size; ++i) {
    int k = i;
    for (int j = k + 1; j < size; ++j)
      if (fabs(temp_matrix.matrix_[k][i]) < fabs(temp_matrix.matrix_[j][i]))
        k = j;

    std::swap(temp_matrix.matrix_[i], temp_matrix.matrix_[k]);

    if (i != k) result = -result;
    result *= temp_matrix.matrix_[i][i];

    for (int j = i + 1; j < size; ++j) {
      temp_matrix.matrix_[i][j] /= temp_matrix.matrix_[i][i];
    }

    for (int j = 0; j < size; ++j) {
      if (j != i && fabs(temp_matrix.matrix_[j][i]) >=
                        std::numeric_limits<double>::epsilon())
        for (k = i + 1; k < size; k++)
          temp_matrix.matrix_[j][k] -=
              temp_matrix.matrix_[j][i] * temp_matrix.matrix_[i][k];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) throw std::logic_error("The matrix is not a square!");

  S21Matrix complementsMatrix = S21Matrix(rows_, cols_);

  if (rows_ == 1 && cols_ == 1) {
    complementsMatrix.matrix_[0][0] = 1;
    return complementsMatrix;
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      complementsMatrix.matrix_[i][j] = CalcComplement(i, j);
    }
  }
  return complementsMatrix;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ != cols_) throw std::logic_error("The matrix is not a square!");
  double determinant = Determinant();
  if (fabs(determinant) < std::numeric_limits<double>::epsilon())
    throw std::runtime_error("The determinant is zero!");
  S21Matrix adjacent_matrix = CalcComplements().Transpose();
  adjacent_matrix.MulNumber(1 / determinant);
  return adjacent_matrix;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &b) {
  if (this == &b) return *this;
  ClearMatrix();
  rows_ = b.rows_;
  cols_ = b.cols_;
  InitMatrix();
  for (int i = 0; i < rows_; ++i) {
    std::copy(b.matrix_[i], b.matrix_[i] + b.cols_, matrix_[i]);
  }
  return *this;
}

void S21Matrix::ClearMatrix() {
  for (int i = 0; i < rows_; ++i) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}

void S21Matrix::InitMatrix() {
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_];
  }
}

//      OPERATORS

double &S21Matrix::operator()(int i, int j) {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    throw std::out_of_range("The index is out of range!");
  return matrix_[i][j];
}

double S21Matrix::operator()(int i, int j) const {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    throw std::out_of_range("The index is out of range!");
  return matrix_[i][j];
}

std::ostream &operator<<(std::ostream &os, const S21Matrix &d) {
  for (int i = 0; i < d.GetRows(); ++i) {
    for (int j = 0; j < d.GetCols(); ++j) {
      std::cout << d(i, j) << '\t';
    }
    std::cout << std::endl;
  }
  return os;
}

S21Matrix operator*(const double &number, const S21Matrix &m1) {
  S21Matrix result = S21Matrix(m1);
  result.MulNumber(number);
  return result;
}

S21Matrix operator+(const S21Matrix &m1, const S21Matrix &m2) {
  S21Matrix result = S21Matrix(m1);
  result.SumMatrix(m2);
  return result;
}
S21Matrix operator-(const S21Matrix &m1, const S21Matrix &m2) {
  S21Matrix result = S21Matrix(m1);
  result.SubMatrix(m2);
  return result;
}
S21Matrix operator*(const S21Matrix &m1, const S21Matrix &m2) {
  S21Matrix result = S21Matrix(m1);
  result.MulMatrix(m2);
  return result;
}

S21Matrix operator*(const S21Matrix &m1, const double &number) {
  S21Matrix result = S21Matrix(m1);
  result.MulNumber(number);
  return result;
}

//      PRIVATE METHODS

double S21Matrix::CalcComplement(int i, int j) {
  S21Matrix minorMatrix = S21Matrix(rows_ - 1, cols_ - 1);
  for (int k = 0; k < rows_; ++k) {
    for (int l = 0; l < cols_; ++l) {
      if (k != i && l != j) {
        int new_i = k > i ? k - 1 : k;
        int new_j = l > j ? l - 1 : l;
        minorMatrix.matrix_[new_i][new_j] = matrix_[k][l];
      }
    }
  }
  double result = minorMatrix.Determinant() * pow(-1, i + j);
  return result;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &b) {
  SumMatrix(b);
  return *this;
}
S21Matrix &S21Matrix::operator-=(const S21Matrix &b) {
  SubMatrix(b);
  return *this;
}
S21Matrix &S21Matrix::operator*=(const S21Matrix &b) {
  MulMatrix(b);
  return *this;
}
S21Matrix &S21Matrix::operator*=(const double &b) {
  MulNumber(b);
  return *this;
}
