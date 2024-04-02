#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() { ResetMatrix(); }

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows > 0 && cols > 0) {
    CreateMatrix(rows, cols);
  } else {
    throw std::invalid_argument("Cols and rows must be bigger than zero");
  }
}

S21Matrix::S21Matrix(const S21Matrix& other) {
  CreateMatrix(other.rows_, other.cols_);
  CopyMatrix(other.matrix_);
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept {
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.ResetMatrix();
}

S21Matrix::~S21Matrix() {
  DeleteMatrix();
  ResetMatrix();
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool res = true;
  if (IsCorrectMatrixesSizes(other) && rows_ == other.rows_ &&
      cols_ == other.cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-6) {
          res = false;
        }
      }
    }
  } else {
    res = false;
  }
  return res;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (IsCorrectMatrixesSizes(other)) {
    if (rows_ == other.rows_ && cols_ == other.cols_) {
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          matrix_[i][j] += other.matrix_[i][j];
        }
      }
    } else {
      throw std::length_error("Lengths of matrixes is not equal");
    }
  } else {
    throw std::invalid_argument("Cols and rows must be bigger than zero");
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (IsCorrectMatrixesSizes(other)) {
    if (rows_ == other.rows_ && cols_ == other.cols_) {
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          matrix_[i][j] -= other.matrix_[i][j];
        }
      }
    } else {
      throw std::length_error("Lengths of matrixes is not equal");
    }
  } else {
    throw std::invalid_argument("Cols and rows must be bigger than zero");
  }
}

void S21Matrix::MulNumber(const double num) {
  if (IsCorrectMatrix_()) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] *= num;
      }
    }
  } else {
    throw std::invalid_argument("Cols and rows must be bigger than zero");
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (IsCorrectMatrixesSizes(other)) {
    if (rows_ == other.cols_ && cols_ == other.rows_) {
      S21Matrix result(rows_, other.cols_);
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < other.cols_; j++) {
          for (int k = 0; k < cols_; ++k) {
            result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
          }
        }
      }
      *this = result;
    } else {
      throw std::length_error(
          "The number of columns of the first matrix is not equal to the "
          "number of rows of the second matrix");
    }
  } else {
    throw std::invalid_argument("Cols and rows must be bigger than zero");
  }
}

S21Matrix S21Matrix::Transpose() {
  if (IsCorrectMatrix_()) {
    S21Matrix Transposed(cols_, rows_);

    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        Transposed.matrix_[j][i] = matrix_[i][j];
      }
    }
    return Transposed;
  } else {
    throw std::invalid_argument("Cols and rows must be bigger than zero");
  }
}

S21Matrix S21Matrix::CalcComplements() {
  if (IsCorrectMatrix_()) {
    if (rows_ == cols_) {
      S21Matrix result(rows_, cols_);
      S21Matrix minor(rows_ - 1, cols_ - 1);

      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          minor = CreateMinor(i, j);
          double minor_det = minor.Determinant();

          double complement = minor_det * ((i + j) % 2 == 0 ? 1 : -1);
          result(i, j) = complement;
        }
      }
      return result;
    } else {
      throw std::invalid_argument("Matrix is now square");
    }
  } else {
    throw std::invalid_argument("Cols and rows must be bigger than zero");
  }
}

double S21Matrix::Determinant() {
  double det = 0.0;

  if (IsCorrectMatrix_()) {
    if (rows_ == cols_) {
      if (rows_ == 2) {
        det = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
      } else if (rows_ == 1) {
        det = matrix_[0][0];
      } else {
        int N = rows_;
        int sign = 1;
        for (int i = 0; i < N; i++) {
          S21Matrix minor(N - 1, N - 1);

          for (int j = 0; j < N - 1; j++) {
            for (int k = 0; k < N - 1; k++) {
              minor.matrix_[j][k] = matrix_[j + 1][(k >= i) ? k + 1 : k];
            }
          }
          double minor_det = minor.Determinant();

          det += sign * matrix_[0][i] * minor_det;
          sign = -sign;
        }
      }
    } else {
      throw std::invalid_argument("Matrix is now square");
    }
  } else {
    throw std::invalid_argument("Cols and rows must be bigger than zero");
  }
  return det;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = 0.0;

  if (IsCorrectMatrix_()) {
    if (rows_ == cols_) {
      det = Determinant();
      if (fabs(det) < 1e-7) {
        throw std::runtime_error("Matrix is singular, inverse does not exist");
      } else {
        S21Matrix result(rows_, cols_);
        if (rows_ == 1) {
          result(0, 0) = 1 / matrix_[0][0];
        } else {
          S21Matrix Complements = CalcComplements();
          S21Matrix transposed_complements = Complements.Transpose();

          for (int i = 0; i < transposed_complements.get_rows(); i++) {
            for (int j = 0; j < transposed_complements.get_cols(); j++) {
              result(i, j) = transposed_complements(i, j) / det;
            }
          }
        }
        return result;
      }
    } else {
      throw std::invalid_argument("Matrix is now square");
    }
  } else {
    throw std::invalid_argument("Cols and rows must be bigger than zero");
  }
}

// overload operators
S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    DeleteMatrix();
    CreateMatrix(other.rows_, other.cols_);
    CopyMatrix(other.matrix_);
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
  if (this != &other) {
    DeleteMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;

    other.ResetMatrix();
  }
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }

void S21Matrix::operator+=(const S21Matrix& other) { SumMatrix(other); }

void S21Matrix::operator-=(const S21Matrix& other) { SubMatrix(other); }

void S21Matrix::operator*=(const S21Matrix& other) { MulMatrix(other); }

void S21Matrix::operator*=(const double num) { MulNumber(num); }

double& S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_ || rows_ < 0 || cols_ < 0) {
    throw std::out_of_range("i and j elements of matrix are out of range");
  }
  return matrix_[i][j];
}

const double& S21Matrix::operator()(int i, int j) const {
  if (i >= rows_ || j >= cols_ || rows_ < 0 || cols_ < 0) {
    throw std::out_of_range("i and j elements of matrix are out     of range");
  }
  return matrix_[i][j];
}

// Another functions
bool S21Matrix::IsCorrectMatrix_() {
  bool result = false;
  if (rows_ > 0 && cols_ > 0 && matrix_ != nullptr) {
    result = true;
  }

  return result;
}

bool S21Matrix::IsCorrectMatrixesSizes(const S21Matrix& other) {
  bool result = false;
  if (IsCorrectMatrix_() &&
      (other.rows_ > 0 && other.cols_ > 0 && other.matrix_ != nullptr)) {
    result = true;
  }

  return result;
}

void S21Matrix::ResetMatrix() {
  cols_ = 0;
  rows_ = 0;
  matrix_ = nullptr;
}

void S21Matrix::CreateMatrix(int rows, int cols) {
  rows_ = rows;
  cols_ = cols;
  matrix_ = new double* [rows_] {};
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]{};
  }
}

void S21Matrix::CopyMatrix(double** matrix) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix[i][j];
    }
  }
}

void S21Matrix::DeleteMatrix() {
  for (int i = 0; i < rows_; ++i) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}

S21Matrix S21Matrix::CreateMinor(int exc_row, int exc_col) {
  int minor_i = 0;
  S21Matrix result(rows_ - 1, cols_ - 1);
  for (int row = 0; row < rows_; row++) {
    if (row == exc_row) {
      continue;
    }

    int minor_j = 0;
    for (int col = 0; col < cols_; col++) {
      if (col == exc_col) {
        continue;
      }

      result.matrix_[minor_i][minor_j] = matrix_[row][col];
      minor_j++;
    }
    minor_i++;
  }
  return result;
}
