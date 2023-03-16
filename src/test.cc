#include <fstream>
#include <gtest/gtest.h>
#include "s21_matrix_oop.h"

static void ReadFromFile(std::ifstream& fstream, S21Matrix *matrix);

class InverseMatrixTests : public ::testing::Test {
 protected:
  S21Matrix identity_matrix { 3, 3 };
  S21Matrix identity_matrix_1x1 { 1, 1 };
  void SetUp() override {
    for (int i { 0 }; i < identity_matrix.get_rows(); ++i) {
      identity_matrix(i, i) = 1;
    }
    identity_matrix_1x1(0, 0) = 1;
  }
};

TEST(ArithmeticTests, SumTest) {
  S21Matrix a { 2, 2 };
  std::ifstream fstr_a { "samples_for_tests/arithmetic_tests/sum2x2.txt" };
  ReadFromFile(fstr_a, &a);
  S21Matrix b { a };
  a += b;
  std::ifstream fstr_b { "samples_for_tests/arithmetic_tests/sum2x2_result.txt" };
  ReadFromFile(fstr_b, &b);
  ASSERT_TRUE(a == b);
}

TEST(ArithmeticTests, SubTest) {
  S21Matrix a { 2, 2 };
  S21Matrix b { a };
  std::ifstream fstr_a { "samples_for_tests/arithmetic_tests/sum2x2.txt" };
  ReadFromFile(fstr_a, &a);
  a -= a;
  ASSERT_TRUE(a == b);
}

TEST(ArithmeticTests, MulNumberTest) {
  S21Matrix a { 2, 2 };
  std::ifstream fstr_a { "samples_for_tests/arithmetic_tests/sum2x2.txt" };
  ReadFromFile(fstr_a, &a);
  S21Matrix b { a };
  a += a;
  b *= 2;
  ASSERT_TRUE(a == b);
}

TEST(ArithmeticTests, MulMatrixTest) {
  S21Matrix a { 3, 2 };
  std::ifstream fstr_a { "samples_for_tests/arithmetic_tests/mult_matrix_3x2.txt" };
  ReadFromFile(fstr_a, &a);
  S21Matrix b { 2, 3 };
  std::ifstream fstr_b { "samples_for_tests/arithmetic_tests/mult_matrix_2x3.txt" };
  ReadFromFile(fstr_b, &b);
  S21Matrix result { 3, 3 };
  std::ifstream fstr_result { "samples_for_tests/arithmetic_tests/mult_matrix_3x3_result.txt" };
  ReadFromFile(fstr_result, &result);
  a *= b;
  ASSERT_TRUE(a == result);
}

TEST(DeterminantTest, SimpleDeterminant) {
  S21Matrix determ { 3, 3 };
  std::ifstream fstr { "samples_for_tests/determ_3x3.txt" };
  ReadFromFile(fstr, &determ);
  ASSERT_DOUBLE_EQ(determ.Determinant(), -107.5);
}

TEST(DeterminantTest, inverted_determ_3x3) {
  S21Matrix determ { 3, 3 };
  std::ifstream fstr { "samples_for_tests/inverted_determ_3x3.txt" };
  ReadFromFile(fstr, &determ);
  ASSERT_DOUBLE_EQ(determ.Determinant(), -18);
}

TEST(DeterminantTest, OneSizeMatrix) {
  S21Matrix one_size_matrix { 1, 1 };
  one_size_matrix(0, 0) = 5;
  ASSERT_DOUBLE_EQ(one_size_matrix.Determinant(), 5);
}

TEST(DeterminantTest, ZeroMatrix) {
  S21Matrix one_size_matrix { 5, 5 };
  ASSERT_DOUBLE_EQ(one_size_matrix.Determinant(), 0);
}

TEST(ResizeMatrix, TestResizing) {
  S21Matrix resizable { 5, 5 };
  resizable.set_rows(1);
  resizable.set_cols(1);
  resizable(0, 0) = 25;
  resizable.set_rows(3);
  resizable.set_cols(3);
  S21Matrix check_matrix { 3, 3 };
  check_matrix(0, 0) = 25;
  ASSERT_TRUE(resizable == check_matrix);
}

TEST(Transpose, TransposeTest) {
  S21Matrix a { 1, 5 };
  a(0, 0) = 5;
  S21Matrix transposed = a.Transpose();
  S21Matrix b { 5, 1 };
  b(0, 0) = 5;
  ASSERT_TRUE(transposed == b);
}

TEST(TestOutput, OutputOverloading) {
  S21Matrix a { 1, 1 };
  std::cout << a;
}

TEST(ConstructorsTests, AssignmentTest) {
  S21Matrix a { 1, 1 };
  a = a;
  ASSERT_TRUE(a == a);
}

TEST(ConstructorsTests, DeepCopyAssignmentTest) {
  S21Matrix a { 1, 1 };
  S21Matrix b { 3, 3 };
  EXPECT_NO_THROW(a = b);
}

TEST(EqualityTest, NotEqualityMatrix) {
  S21Matrix a { 5, 1 };
  S21Matrix b { a };
  b(3, 0) = 3.4;
  ASSERT_FALSE(a == b);
}

TEST(EqualityTest, DiffentSize) {
  S21Matrix a { 5, 1 };
  S21Matrix b { 3, 3 };
  ASSERT_FALSE(a == b);
}

TEST(MatrixComplementTests, Matrix3x3Complements) {
  S21Matrix test { 3, 3 };
  std::ifstream fstr_test { "samples_for_tests/calc_complements/3x3.txt" };
  ReadFromFile(fstr_test, &test);
  S21Matrix result { 3, 3 };
  std::ifstream fstr_result { "samples_for_tests/calc_complements/3x3_result.txt" };
  ReadFromFile(fstr_result, &result);
  ASSERT_TRUE(test.CalcCompements() == result);
}

TEST(MatrixComplementTests, Matrix1x1Complements) {
  S21Matrix test { 1, 1 };
  test(0, 0) = 45;
  S21Matrix result { 1, 1 };
  result(0, 0) = 45;
  ASSERT_TRUE(test.CalcCompements() == result);
}

TEST(MinorTests, Minor3x3) {
  S21Matrix a { 3, 3 };
  S21Matrix result { 2, 2 };
  ASSERT_TRUE(result == a.Minor(0, 0));
}

TEST(MinorTests, Minor1x1) {
  S21Matrix one { 1, 1 };
  ASSERT_TRUE(one == one.Minor(0, 0));
}

TEST(MinorTests, NegativeRangeError) {
  S21Matrix one { 1, 1 };
  EXPECT_THROW(one == one.Minor(-1, 0), std::range_error);
  EXPECT_THROW(one == one.Minor(0, -1), std::range_error);
}

TEST_F(InverseMatrixTests, Inverse3x3) {
  S21Matrix test { 3, 3 };
  std::ifstream fstr_test { "samples_for_tests/determ_3x3.txt" };
  ReadFromFile(fstr_test, &test);
  S21Matrix inverse { test.InverseMatrix() };
  test.MulMatrix(inverse);
  ASSERT_TRUE(test == identity_matrix);
}

TEST_F(InverseMatrixTests, Inverse1x1) {
  S21Matrix test { 1, 1 };
  test(0, 0) = 45;
  S21Matrix inverse { test.InverseMatrix() };
  test.MulMatrix(inverse);
  ASSERT_TRUE(test == identity_matrix_1x1);
}

TEST(BinaryOverloading, SumMatrix) {
  S21Matrix a { 3, 3 };
  a(1, 1) = 5;
  S21Matrix b { 3, 3 };
  b(1, 1) = 10;
  S21Matrix result{ 3, 3 };
  result(1, 1) = 15;
  ASSERT_TRUE(a + b == result);
}

TEST(BinaryOverloading, SubMatrix) {
  S21Matrix a { 3, 3 };
  a(1, 1) = 5;
  S21Matrix b { 3, 3 };
  b(1, 1) = 10;
  S21Matrix result{ 3, 3 };
  result(1, 1) = -5;
  ASSERT_TRUE(a - b == result);
}

TEST(BinaryOverloading, MulMatrixByNumber) {
  S21Matrix a { 3, 3 };
  a(1, 1) = 5;
  S21Matrix result{ 3, 3 };
  result(1, 1) = 15;
  ASSERT_TRUE(a * 3 == result);
  ASSERT_TRUE(3 * a == result);
}

TEST(BinaryOverloading, MulMatrixTest) {
  S21Matrix a { 3, 2 };
  std::ifstream fstr_a { "samples_for_tests/arithmetic_tests/mult_matrix_3x2.txt" };
  ReadFromFile(fstr_a, &a);
  S21Matrix b { 2, 3 };
  std::ifstream fstr_b { "samples_for_tests/arithmetic_tests/mult_matrix_2x3.txt" };
  ReadFromFile(fstr_b, &b);
  S21Matrix result { 3, 3 };
  std::ifstream fstr_result { "samples_for_tests/arithmetic_tests/mult_matrix_3x3_result.txt" };
  ReadFromFile(fstr_result, &result);
  ASSERT_TRUE(a * b == result);
}

TEST(Exceptions, NonSquareDeterminant) {
  S21Matrix a { 1, 3 };
  ASSERT_THROW(a.Determinant(), std::logic_error);
}

TEST(Exceptions, NegativeRangeMatrices) {
  EXPECT_THROW((S21Matrix{ -1, 3 }), std::range_error);
  EXPECT_THROW((S21Matrix{ 1, -3 }), std::range_error);
}

TEST(Exceptions, NegativeResizing) {
  S21Matrix a {};
  EXPECT_THROW(a.set_rows(-5), std::range_error);
  EXPECT_THROW(a.set_cols(-5), std::range_error);
}

TEST(Exceptions, OutOfRangeAccess) {
  S21Matrix a { 3, 3 };
  EXPECT_THROW(a(-5, 0), std::range_error);
  S21Matrix b { 3, 3 };
  EXPECT_THROW(a(0, -5), std::range_error);
}

TEST(Exceptions, ConstOutOfRangeAccess) {
  const S21Matrix a { 3, 3 };
  EXPECT_THROW(a(-5, 0), std::range_error);
  const S21Matrix b { 3, 3 };
  EXPECT_THROW(a(0, -5), std::range_error);
}

TEST(Exceptions, SumDifferentSize) {
  EXPECT_THROW((S21Matrix{ 3, 3 } + S21Matrix{ 2, 5 }), std::logic_error);
}

TEST(Exceptions, SubDifferentSize) {
  EXPECT_THROW((S21Matrix{ 3, 3 }.SubMatrix(S21Matrix{ 2, 5 })), std::logic_error);
}

TEST(Exceptions, MulDifferentSize) {
  EXPECT_THROW((S21Matrix{ 3, 3 } * S21Matrix{ 2, 5 }), std::logic_error);
}

TEST(Exceptions, NonSquareComplements) {
  EXPECT_THROW((S21Matrix{ 2, 3 }.CalcCompements()), std::logic_error);
}

TEST(Exceptions, NonSquareMinor) {
  EXPECT_THROW((S21Matrix{ 2, 3 }.Minor(0, 0)), std::logic_error);
}

TEST(Exceptions, NonSquareInverseMatrix) {
  EXPECT_THROW((S21Matrix{ 2, 3 }.InverseMatrix()), std::logic_error);
}

TEST(Exceptions, ZeroDeterminantInverseMatrix) {
  EXPECT_THROW((S21Matrix{ 3, 3 }.InverseMatrix()), std::logic_error);
}

TEST(MoveTest, MoveAssignment) {
  S21Matrix a { 3, 3 };
  EXPECT_NO_THROW(a = std::move(a));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

static void ReadFromFile(std::ifstream& fstream, S21Matrix *matrix) {
  for (int i { 0 }; i < matrix->get_rows(); ++i) {
    for (int j { 0 }; j < matrix->get_cols(); ++j) {
      fstream >> (*matrix)(i, j);
    }
  }
}

