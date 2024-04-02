#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <cmath>
#include <stdexcept>
class S21Matrix {
 private:
  // Attributes
  int rows_, cols_;  // Rows and columns
  double** matrix_;  // Pointer to the memory where the matrix is allocated

  void ResetMatrix();

  void CreateMatrix(int rows, int cols);

  void CopyMatrix(double** matrix);
  
  void DeleteMatrix();

  bool IsCorrectMatrix_();

  bool IsCorrectMatrixesSizes(const S21Matrix& other);

  S21Matrix CreateMinor(int exc_row, int exc_col);

 public:
  S21Matrix();  // Default constructor

  S21Matrix(int rows, int cols);

  S21Matrix(const S21Matrix& other);

  S21Matrix(S21Matrix&& other) noexcept;

  ~S21Matrix();  // Destructor

  int get_cols() const { return cols_; };

  int get_rows() const { return rows_; };

  bool EqMatrix(const S21Matrix& other);

  void SumMatrix(const S21Matrix& other);

  void SubMatrix(const S21Matrix& other);

  void MulNumber(const double num);

  void MulMatrix(const S21Matrix& other);

  S21Matrix Transpose();

  S21Matrix CalcComplements();

  double Determinant();

  S21Matrix InverseMatrix();

  // overload operators
  S21Matrix& operator=(const S21Matrix& other);

  S21Matrix& operator=(S21Matrix&& other) noexcept;

  S21Matrix operator+(const S21Matrix& other);

  S21Matrix operator-(const S21Matrix& other);

  S21Matrix operator*(const S21Matrix& other);

  S21Matrix operator*(const double num);

  bool operator==(const S21Matrix& other);

  void operator+=(const S21Matrix& other);

  void operator-=(const S21Matrix& other);

  void operator*=(const S21Matrix& other);

  void operator*=(const double num);

  double& operator()(int i, int j);
  
  const double& operator()(int i, int j) const;
};

#endif