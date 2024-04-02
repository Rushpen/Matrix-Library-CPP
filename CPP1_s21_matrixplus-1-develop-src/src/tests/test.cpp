#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

class MatrixTest : public ::testing::Test {
 public:
  S21Matrix FillWithStep(int rows, int cols, double elem) {
    S21Matrix tmp(rows, cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        tmp(i, j) = elem;
        elem++;
      }
    }
    return tmp;
  }

  S21Matrix FillMatrix(int rows, int cols, double elem) {
    S21Matrix tmp(rows, cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        tmp(i, j) = elem;
      }
    }
    return tmp;
  }

  S21Matrix *matrix_a;
  S21Matrix *matrix_b;
  S21Matrix *matrix_c;
  S21Matrix *matrix_d;
  S21Matrix *matrix_e;
  S21Matrix *matrix_f;
  S21Matrix *matrix_g;
  S21Matrix *matrix_h;
  S21Matrix *matrix_b_plus_h;
  S21Matrix *matrix_i;

  void SetUp() override {
    matrix_a = new S21Matrix(3, 2);
    matrix_b = new S21Matrix(2, 3);
    *matrix_b = FillWithStep(2, 3, 1.1234567);
    matrix_c = new S21Matrix(2, 3);
    *matrix_c = FillWithStep(2, 3, 1.1234567);
    matrix_d = new S21Matrix(2, 3);
    *matrix_d = FillWithStep(2, 3, 56.987654321);
    matrix_e = new S21Matrix(3, 3);
    matrix_f = new S21Matrix(1, 1);
    matrix_g = new S21Matrix(4, 2);
    matrix_h = new S21Matrix(2, 3);
    *matrix_h = FillMatrix(2, 3, 1);
    matrix_b_plus_h = new S21Matrix(2, 3);
    *matrix_b_plus_h = FillWithStep(2, 3, 2.1234567);
    matrix_i = new S21Matrix(3, 2);
    *matrix_i = FillWithStep(3, 2, 56.987654321);
  }

  void TearDown() override {
    delete matrix_a;
    delete matrix_b;
    delete matrix_c;
    delete matrix_d;
    delete matrix_e;
    delete matrix_f;
    delete matrix_g;
    delete matrix_h;
    delete matrix_b_plus_h;
    delete matrix_i;
  }
};

TEST_F(MatrixTest, Testget_rows) { ASSERT_EQ(3, matrix_a->get_rows()); }

TEST_F(MatrixTest, Testget_cols) { ASSERT_EQ(2, matrix_a->get_cols()); }

TEST_F(MatrixTest, TestEqual) { ASSERT_TRUE(matrix_b->EqMatrix(*matrix_c)); }

TEST_F(MatrixTest, TestEqualDifferentValues) {
  ASSERT_FALSE(matrix_b->EqMatrix(*matrix_d));
}

TEST_F(MatrixTest, TestEqualDifferentSize) {
  ASSERT_FALSE(matrix_b->EqMatrix(*matrix_g));
}

TEST_F(MatrixTest, TestEqualDifferentValues2) {
  (*matrix_c)(1, 2) = 1.1234568;
  ASSERT_FALSE(matrix_b->EqMatrix(*matrix_c));
}

TEST_F(MatrixTest, TestSum) {
  matrix_b->SumMatrix(*matrix_h);
  ASSERT_TRUE(matrix_b->EqMatrix(*matrix_b_plus_h));
}

TEST_F(MatrixTest, TestSub) {
  matrix_b_plus_h->SubMatrix(*matrix_h);
  ASSERT_TRUE(matrix_b_plus_h->EqMatrix(*matrix_b));
}

TEST_F(MatrixTest, TestMultMatrix) {
  matrix_c->MulMatrix(*matrix_i);
  S21Matrix matrix_correct(2, 2);
  matrix_correct(0, 0) = 379.773189;
  matrix_correct(0, 1) = 386.143559;
  matrix_correct(1, 0) = 910.662078;
  matrix_correct(1, 1) = 926.032448;
  ASSERT_TRUE(matrix_c->EqMatrix(matrix_correct));
}

TEST_F(MatrixTest, TestMultNum) {
  matrix_c->MulNumber(1.12345670);
  S21Matrix matrix_correct(2, 3);
  matrix_correct(0, 0) = 1.262155;
  matrix_correct(0, 1) = 2.385612;
  matrix_correct(0, 2) = 3.509068;
  matrix_correct(1, 0) = 4.632525;
  matrix_correct(1, 1) = 5.755982;
  matrix_correct(1, 2) = 6.879438;
  ASSERT_TRUE(matrix_c->EqMatrix(matrix_correct));
}

TEST_F(MatrixTest, TestTranspose) {
  S21Matrix matrix_correct(3, 2);
  matrix_correct(0, 0) = 1.123457;
  matrix_correct(1, 0) = 2.123457;
  matrix_correct(2, 0) = 3.123457;
  matrix_correct(0, 1) = 4.123457;
  matrix_correct(1, 1) = 5.123457;
  matrix_correct(2, 1) = 6.123457;
  ASSERT_TRUE(matrix_c->Transpose().EqMatrix(matrix_correct));
}

TEST_F(MatrixTest, TestCalcComp) {
  (*matrix_e)(0, 0) = 1;
  (*matrix_e)(0, 1) = 2;
  (*matrix_e)(0, 2) = 3;
  (*matrix_e)(1, 0) = 0;
  (*matrix_e)(1, 1) = 4;
  (*matrix_e)(1, 2) = 2;
  (*matrix_e)(2, 0) = 5;
  (*matrix_e)(2, 1) = 2;
  (*matrix_e)(2, 2) = 1;
  S21Matrix matrix_correct(3, 3);
  matrix_correct(0, 0) = 0;
  matrix_correct(0, 1) = 10;
  matrix_correct(0, 2) = -20;
  matrix_correct(1, 0) = 4;
  matrix_correct(1, 1) = -14;
  matrix_correct(1, 2) = 8;
  matrix_correct(2, 0) = -8;
  matrix_correct(2, 1) = -2;
  matrix_correct(2, 2) = 4;
  ASSERT_TRUE(matrix_e->CalcComplements().EqMatrix(matrix_correct));
}

TEST_F(MatrixTest, TestDeterminant) {
  (*matrix_e)(0, 0) = 1;
  (*matrix_e)(0, 1) = 2;
  (*matrix_e)(0, 2) = -4;
  (*matrix_e)(1, 0) = 2;
  (*matrix_e)(1, 1) = 3;
  (*matrix_e)(1, 2) = 1;
  (*matrix_e)(2, 0) = 3;
  (*matrix_e)(2, 1) = -1;
  (*matrix_e)(2, 2) = 2;
  ASSERT_DOUBLE_EQ(49, matrix_e->Determinant());
}

TEST_F(MatrixTest, TestDeterminant2) {
  (*matrix_f)(0, 0) = 1.56;
  ASSERT_DOUBLE_EQ(1.56, matrix_f->Determinant());
}

TEST_F(MatrixTest, TestInverse) {
  (*matrix_e)(0, 0) = 1;
  (*matrix_e)(0, 1) = -2;
  (*matrix_e)(0, 2) = 3;
  (*matrix_e)(1, 0) = 0;
  (*matrix_e)(1, 1) = 4;
  (*matrix_e)(1, 2) = -1;
  (*matrix_e)(2, 0) = 5;
  (*matrix_e)(2, 1) = 0;
  (*matrix_e)(2, 2) = 0;
  S21Matrix matrix_correct(3, 3);
  matrix_correct(0, 0) = 0;
  matrix_correct(0, 1) = 0;
  matrix_correct(0, 2) = 0.2;
  matrix_correct(1, 0) = 0.1;
  matrix_correct(1, 1) = 0.3;
  matrix_correct(1, 2) = -0.02;
  matrix_correct(2, 0) = 0.4;
  matrix_correct(2, 1) = 0.2;
  matrix_correct(2, 2) = -0.08;
  ASSERT_TRUE(matrix_e->InverseMatrix().EqMatrix(matrix_correct));
}

TEST_F(MatrixTest, TestInverse2) {
  (*matrix_f)(0, 0) = 1.56;
  ASSERT_DOUBLE_EQ(1.0 / 1.56, matrix_f->InverseMatrix()(0, 0));
}

TEST_F(MatrixTest, TestSumOperator) {
  S21Matrix matrix_sum = *matrix_b + *matrix_h;
  ASSERT_TRUE(matrix_sum.EqMatrix(*matrix_b_plus_h));
}

TEST_F(MatrixTest, TestSubOperator) {
  S21Matrix matrix_sub = *matrix_b_plus_h - *matrix_h;
  ASSERT_TRUE(matrix_sub.EqMatrix(*matrix_b));
}

TEST_F(MatrixTest, TestMultMatrixOperator) {
  S21Matrix matrix_mult = (*matrix_b) * (*matrix_i);
  S21Matrix matrix_correct(2, 2);
  matrix_correct(0, 0) = 379.773189;
  matrix_correct(0, 1) = 386.143559;
  matrix_correct(1, 0) = 910.662078;
  matrix_correct(1, 1) = 926.032448;
  ASSERT_TRUE(matrix_mult.EqMatrix(matrix_correct));
}

TEST_F(MatrixTest, TestMultNumOperator) {
  S21Matrix matrix_mult = (*matrix_b) * 1.12345670;
  S21Matrix matrix_correct(2, 3);
  matrix_correct(0, 0) = 1.262155;
  matrix_correct(0, 1) = 2.385612;
  matrix_correct(0, 2) = 3.509068;
  matrix_correct(1, 0) = 4.632525;
  matrix_correct(1, 1) = 5.755982;
  matrix_correct(1, 2) = 6.879438;
  ASSERT_TRUE(matrix_mult.EqMatrix(matrix_correct));
}

TEST_F(MatrixTest, TestEqualOperator) {
  ASSERT_TRUE((*matrix_b) == (*matrix_c));
}

TEST_F(MatrixTest, TestSumEqOperator) {
  (*matrix_b) += (*matrix_h);
  ASSERT_TRUE(matrix_b->EqMatrix(*matrix_b_plus_h));
}

TEST_F(MatrixTest, TestSubEqOperator) {
  (*matrix_b_plus_h) -= (*matrix_h);
  ASSERT_TRUE(matrix_b_plus_h->EqMatrix(*matrix_b));
}

TEST_F(MatrixTest, TestMultEqMatrixOperator) {
  (*matrix_b) *= (*matrix_i);
  S21Matrix matrix_correct(2, 2);
  matrix_correct(0, 0) = 379.773189;
  matrix_correct(0, 1) = 386.143559;
  matrix_correct(1, 0) = 910.662078;
  matrix_correct(1, 1) = 926.032448;
  ASSERT_TRUE(matrix_b->EqMatrix(matrix_correct));
}

TEST_F(MatrixTest, TestMultEqNumOperator) {
  (*matrix_b) *= 1.12345670;
  S21Matrix matrix_correct(2, 3);
  matrix_correct(0, 0) = 1.262155;
  matrix_correct(0, 1) = 2.385612;
  matrix_correct(0, 2) = 3.509068;
  matrix_correct(1, 0) = 4.632525;
  matrix_correct(1, 1) = 5.755982;
  matrix_correct(1, 2) = 6.879438;
  ASSERT_TRUE(matrix_b->EqMatrix(matrix_correct));
}

TEST_F(MatrixTest, TestMoveAssignment) {
  S21Matrix matrix_moved = std::move(*matrix_b);
  ASSERT_TRUE(matrix_moved.EqMatrix(*matrix_c));
}

TEST_F(MatrixTest, TestIndexOperatorConst) {
  const S21Matrix matrix_const(2, 2);
  ASSERT_EQ(0, matrix_const(1, 1));
}

TEST_F(MatrixTest, TestParameterizedConstructorException) {
  ASSERT_ANY_THROW(S21Matrix empty(-1, 2));
}

TEST_F(MatrixTest, TestSumMatrixException) {
  ASSERT_ANY_THROW(matrix_a->SumMatrix(*matrix_f));
}

TEST_F(MatrixTest, TestSubMatrixException) {
  ASSERT_ANY_THROW(matrix_a->SubMatrix(*matrix_f));
}

TEST_F(MatrixTest, TestMulMatrixException) {
  ASSERT_ANY_THROW(matrix_a->MulMatrix(*matrix_f));
}

TEST_F(MatrixTest, TestCalComplementsException) {
  ASSERT_ANY_THROW(matrix_a->CalcComplements());
}

TEST_F(MatrixTest, TestDeterminantException) {
  ASSERT_ANY_THROW(matrix_a->Determinant());
}

TEST_F(MatrixTest, TestInverseMatrixException1) {
  ASSERT_ANY_THROW(matrix_a->InverseMatrix());
}

TEST_F(MatrixTest, TestInverseMatrixException2) {
  ASSERT_ANY_THROW(matrix_f->InverseMatrix());
}

TEST_F(MatrixTest, TestIndexOperatorException) {
  ASSERT_ANY_THROW((*matrix_f)(3, -2));
}

TEST_F(MatrixTest, TestIndexOperatorConstException) {
  const S21Matrix matrix_const(2, 2);
  ASSERT_ANY_THROW(matrix_const(3, -2));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}