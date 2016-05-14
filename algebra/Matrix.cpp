#include "Matrix.hpp"

#include <sstream>

Matrix::Matrix() :
	m_rows(3), 
	m_cols(3), 
	m_values(NULL)
{
	init();
}



Matrix::Matrix(unsigned int rows, unsigned int cols) :
	m_rows(rows),
	m_cols(cols),
	m_values(NULL)
{
	init();
}



Matrix::Matrix(const Matrix& matrix) :
	m_rows(0), 
	m_cols(0),
	m_values(NULL)
{
	copy(matrix);
}

/**
 * destructor
 */
Matrix::~Matrix()
{
  delete[] m_values;
}

/**
 * @param row
 *     the row
 * @param col
 *     the column
 *
 * @return reference of given position
 * @throw Exception, if col or row index is out of bounds
 */
double&
Matrix::operator()(unsigned int row, unsigned int col) throw (Exception)
{
  if (col >= m_cols || row >= m_rows)
    {
      throw Exception(__FILE__, __LINE__, "Error, index out of bounds");
    }
  else
    {
      return *(m_values + row * m_cols + col);
    }
}

/**
 * self assignment safe assignment operator
 *
 * @param matrix
 *     a matrix to assign
 *
 * @return a reference to this
 */
Matrix&
Matrix::operator=(const Matrix &matrix)
{
  copy(matrix);
  return *this;
}

/**
 * @param matrix
 *     a matrix to compare
 *
 * @return true if matrices are equal
 */
bool
Matrix::operator==(const Matrix& matrix) const
{
  if (!matchesDimensions(matrix))
    {
      return false;
    }

  for (unsigned int i = 0; i < (m_rows * m_cols); i++)
    {
      if (m_values[i] != matrix.m_values[i])
        {
          return false;
        }
    }

  return true;
}

/**
 * @param matrix
 *     a matrix to compare
 *
 * @return true if matrices are not equal
 */
bool
Matrix::operator!=(const Matrix& matrix) const
{
  return !(*this == matrix);
}

/**
 * @param matrix
 *     the matrix to add
 *
 * @return sum matrix this+m
 * @throw Exception if matrix sizes does not match
 */
Matrix
Matrix::operator+(const Matrix& matrix) const throw (Exception)
{
  if (m_cols != matrix.m_cols || m_rows != matrix.m_rows)
    {
      throw Exception(__FILE__, __LINE__,
          "Error, trying to add matrices with different dimensions");
    }
  else
    {
      Matrix sum(m_rows, m_cols);
      for (unsigned int row = 0; row < m_rows; row++)
        {
          for (unsigned int col = 0; col < m_cols; col++)
            {
              sum(row, col) = get(row, col) + matrix.get(row, col);
            }
        }
      return sum;
    }
}

/**
 * @param matrix
 *     another matrix to add
 *
 * @return *this
 */
Matrix&
Matrix::operator+=(const Matrix& matrix)
{
  (*this) = (*this) + matrix;
  return (*this);
}

/**
 * @param matrix
 *     a matrix to substract
 *
 * @return the difference this - m
 */
Matrix
Matrix::operator-(const Matrix& matrix) const
{
  return (*this) + (-1 * matrix);
}

/**
 * @param matrix
 *     a matrix
 *
 * @return a matrix where all elements have inverse signs
 */
Matrix
operator-(const Matrix& matrix)
{
  return (-1 * matrix);
}

/**
 * @param matrix
 *     a matrix to substract
 *
 * @return *this
 */
Matrix&
Matrix::operator-=(const Matrix& matrix)
{
  (*this) = (*this) - matrix;
  return (*this);
}

/**
 * @param factor
 *     a factor to multiply
 *
 * @return a matrix where all elements are multiplied with the given factor
 */
Matrix
Matrix::operator*(double factor) const
{
  Matrix product(m_rows, m_cols);

  for (unsigned int row = 0; row < m_rows; row++)
    {
      for (unsigned int col = 0; col < m_cols; col++)
        {
          product(row, col) = factor * get(row, col);
        }
    }

  return product;
}

/**
 * @param factor
 *     a factor to multiply
 * @param matrix
 *     a matrix to multiply with f
 *
 * @return a matrix where all elements are multiplied with the given factor
 */
Matrix
operator*(double factor, const Matrix& matrix)
{
  return matrix * factor;
}

/**
 * @param matrix
 *     a matrix to multiply
 *
 * @return the matrix product this*matrix
 * @throw Exception if matrix dimensions do not match
 */
Matrix
Matrix::operator*(const Matrix& matrix) const throw (Exception)
{
  if (m_cols != matrix.m_rows)
    {
      stringstream ss;
      ss << "Error, matrix multiplication of " << m_rows << "x" << m_cols
          << " and " << matrix.m_rows << "x" << matrix.m_cols
          << " is not defined";
      throw Exception(__FILE__, __LINE__, ss.str());
    }

  Matrix product(m_rows, matrix.m_cols);
  for (unsigned int row = 0; row < product.m_rows; row++)
    {
      for (unsigned int col = 0; col < product.m_cols; col++)
        {
          for (unsigned int i = 0; i < m_cols; i++)
            {
              product(row, col) += get(row, i) * matrix.get(i, col);
            }
        }
    }

  return product;
}

/**
 * @param vector
 *     a vector to multiplicate
 * 
 * @return the result of this*vector
 * @throw Exception if dimensions do not match
 */
Vector
Matrix::operator*(const Vector& vector) const throw (Exception)
{
  if (m_cols != vector.dim())
    {
      throw Exception(__FILE__, __LINE__,
          "Error, vector dimension must be the same as the number of matrix columns");
    }

  Vector product(m_rows);
  for (unsigned int row = 0; row < m_rows; row++)
    {
      for (unsigned int col = 0; col < m_cols; col++)
        {
          product[row] += get(row, col) * vector.m_values.at(col);
        }
    }
  return product;
}

/**
 * @param matrix
 *     another matrix to multiply
 * 
 * @return *this
 * @throw Exceptions if dimensions does not match
 */
Matrix&
Matrix::operator*=(const Matrix& matrix) throw (Exception)
{
  (*this) = (*this) * matrix;
  return (*this);
}

/**
 * @param divisor
 *     the divisor to divide the matrix elements
 *
 * @return this/divisor
 * @throw Exception if divisor=0 (division by zero)
 */
Matrix
Matrix::operator/(double divisor) const throw (Exception)
{
  if (divisor == 0)
    {
      throw Exception(__FILE__, __LINE__, "Error, division by zero");
    }
  return (1 / divisor) * (*this);
}

/**
 * @param matrix
 *     the divisor to divide the matrix 
 *
 * @return this * matrix^-1
 * @throw Exception if matrix is not invertible (division by zero) or the dimensions do not match
 */
Matrix
Matrix::operator/(const Matrix& matrix) const throw (Exception)
{
  if (!matrix.isInvertible())
    {
      throw Exception(__FILE__, __LINE__, "Error, operand is not invertible");
    }
  return (*this) * matrix.inverse();
}

/**
 * @param os
 *     an ofstream object
 * @param matrix
 *     a matrix to stream
 *
 * @return the ofstream object containing the matrix data
 */
ostream&
operator<<(ostream& os, const Matrix& matrix)
{
  if (matrix.m_rows == 0 || matrix.m_cols == 0)
    {
      return os << "Matrix" << matrix.m_rows << "X" << matrix.m_cols
          << "(empty)";
    }
  else
    {
      os << "Matrix" << matrix.m_rows << "X" << matrix.m_cols << "\n";
      for (unsigned int row = 0; row < matrix.m_rows; row++)
        {
          os << " (" << matrix.get(row, 0);
          for (unsigned int col = 1; col < matrix.m_cols; col++)
            {
              os << ", " << matrix.get(row, col);
            }
          os << ")\n";
        }
      return os;
    }
}

/**
 * @return the number of rows
 */
unsigned int
Matrix::getRows() const
{
  return m_rows;
}

/**
 * @return the number of columns
 */
unsigned int
Matrix::getCols() const
{
  return m_cols;
}

/**
 * @param row
 *     the row [0...]
 * @param col
 *     the column [0...]
 *
 * @return the value at the given row/col position
 * @throw Exception if index is out of bounds
 */
double
Matrix::get(unsigned int row, unsigned int col) const throw (Exception)
{
  if (col >= m_cols || row >= m_rows)
    {
      throw Exception(__FILE__, __LINE__, "Error, index out of bounds");
    }
  else
    {
      return m_values[row * m_cols + col];
    }
}

/**
 * @param row
 *     the row [0...]
 * @param col
 *     the column [0...]
 * @param value
 *     the value to set
 *
 * @throw Exception if index is out of bounds
 */
void
Matrix::set(unsigned int row, unsigned int col, double value) throw (Exception)
{
  if (col >= m_cols || row >= m_rows)
    {
      throw Exception(__FILE__, __LINE__, "Error, index out of bounds");
    }
  m_values[row * m_cols + col] = value;
}

/**
 * @param remove_row
 *     the row that should be removed
 * @param remove_col
 *     the column that should be removed
 *
 * @return a submatrix, without the given row and column
 * @throw Exception if row or col is out of bounds
 */
Matrix
Matrix::submatrix(unsigned int remove_row, unsigned int remove_col) const
    throw (Exception)
{
  if (remove_col >= m_cols || remove_row >= m_rows)
    {
      throw Exception(__FILE__, __LINE__, "Error, index out of bounds");
    }

  Matrix submatrix(m_rows - 1, m_cols - 1);

  for (unsigned int row = 0, src_row = 0; row < submatrix.m_rows; row++, src_row++)
    {
      for (unsigned int col = 0, src_col = 0; col < submatrix.m_cols; col++, src_col++)
        {
          if (src_row == remove_row)
            {
              src_row++;
            }
          if (src_col == remove_col)
            {
              src_col++;
            }

          submatrix(row, col) = get(src_row, src_col);
        }
    }

  return submatrix;
}

/**
 * @return the transposed matrix
 */
Matrix
Matrix::transposed() const
{
  Matrix transposed_matrix(m_cols, m_rows);

  for (unsigned int row = 0; row < m_rows; row++)
    {
      for (unsigned int col = 0; col < m_cols; col++)
        {
          transposed_matrix(col, row) = get(row, col);
        }
    }

  return transposed_matrix;
}

/**
 * @param matrix
 *     a matrix to transpose
 *
 * @return the transposed matrix
 */
Matrix
transposed(const Matrix& matrix)
{
  return matrix.transposed();
}

/**
 * @param row
 *     the row index [0...]
 * @param col
 *     the column index [0...]
 *
 * @return the cofactor (-1)^(row+col)* Minor(row,col)
 * @throw Exception if row or column index is out of bounds
 */
double
Matrix::cofactor(unsigned int row, unsigned int col) const throw (Exception)
{
  if (!isSquare())
    {
      throw Exception(__FILE__, __LINE__,
          "Error, there are no cofactors in non square matrices");
    }
  if (col >= m_cols || row >= m_rows)
    {
      throw Exception(__FILE__, __LINE__, "Error, index out of bounds");
    }

  return pow(double(-1), double(row + col)) * submatrix(row, col).det();
}

/**
 * calculates the adjugate matrix
 * transposing of cofactor matrix will be made without calling transposed()
 *
 * @return the adjugate matrix
 * @throw Exception if this is no square matrix
 */
Matrix
Matrix::adjugate() const throw (Exception)
{
  if (!isSquare())
    {
      throw Exception(__FILE__, __LINE__,
          "Error, there are no adjugates of non square matrices");
    }

  Matrix adjugate_matrix(m_rows, m_cols);

  for (unsigned int row = 0; row < m_rows; row++)
    {
      for (unsigned int col = 0; col < m_cols; col++)
        {
          adjugate_matrix(row, col) = cofactor(col, row);
        }
    }

  return adjugate_matrix;
}

/**
 * calculates the adjugate matrix
 * transposing of cofactor matrix will be made without calling transposed()
 *
 * @param matrix
 *     a matrix to calculate the adjugate matrix from
 *
 * @return the adjugate matrix
 * @throw Exception if this is no square matrix
 */
Matrix
adjugate(const Matrix& matrix) throw (Exception)
{
  return matrix.adjugate();
}

/**
 * at the moment, this matrix class operates on real numbers
 * and in this case the conjugate is equal to the transposed
 *
 * @return the conjugate transposed matrix of this
 */
Matrix
Matrix::conjugate() const
{
  return transposed();
}

/**
 * at the moment, this matrix class operates on real numbers
 * and in this case the conjugate is equal to the transposed
 *
 * @param matrix
 *     a matrix to conjugate
 *
 * @return the conjugate transposed matrix
 */
Matrix
conjugate(const Matrix& matrix)
{
  return matrix.conjugate();
}

/**
 * @return the inverse matrix
 * @throw Exception if this is no square matrix or this is not invertible
 */
Matrix
Matrix::inverse() const throw (Exception)
{
  if (!isSquare())
    {
      throw Exception(__FILE__, __LINE__,
          "Error, there are no inverse on non square matrices");
    }
  if (!isInvertible())
    {
      throw Exception(__FILE__, __LINE__,
          "Error, matrix is not invertible (determinant=0)");
    }

  return (1 / det()) * adjugate();
}

/**
 * @param matrix
 *     a matrix to invert
 *
 * @return the inverse matrix
 * @throw Exception if this is no square matrix or this is not invertible
 */
Matrix
inverse(const Matrix& matrix) throw (Exception)
{
  return matrix.inverse();
}

/**
 * @return the determinant of this matrix
 * @throw Exception if this is no square matrix
 */
double
Matrix::det() const throw (Exception)
{
  if (!isSquare())
    {
      stringstream ss;
      ss << "Error, calculating the determinant of a " << m_rows << "x"
          << m_cols << " matrix is not defined";
      throw Exception(__FILE__, __LINE__, ss.str());
    }

  if (m_rows == 1)
    {
      return m_values[0];
    }
  else if (m_rows == 2)
    {
      return get(0, 0) * get(1, 1) - get(0, 1) * get(1, 0);
    }
  else if (m_rows == 3)
    {
      return get(0, 0) * get(1, 1) * get(2, 2) + get(0, 1) * get(1, 2) * get(2,
          0) + get(0, 2) * get(1, 0) * get(2, 1) - get(0, 2) * get(1, 1) * get(
          2, 0) - get(0, 1) * get(1, 0) * get(2, 2) - get(0, 0) * get(1, 2)
          * get(2, 1);
    }
  else
    {
      double determinant = 0;
      for (unsigned int col = 0; col < m_cols; col++)
        {
          determinant += get(0, col) * cofactor(0, col);
        }
      return determinant;
    }
}

/**
 * @return the trace of this matrix
 * @throw Exception if matrix is not square
 */
double
Matrix::trace() const throw (Exception)
{
  if (!isSquare())
    {
      throw Exception(__FILE__, __LINE__,
          "Error, trace is only defined on square matrices");
    }

  double trace = 0;
  for (unsigned int i = 0; i < m_cols; i++)
    {
      trace += get(i, i);
    }
  return trace;
}

/**
 * @param matrix
 *     a matrix to calculate the trace
 *
 * @return the trace of this matrix
 * @throw Exception if matrix is not square
 */
double
trace(const Matrix& matrix) throw (Exception)
{
  return matrix.trace();
}

/**
 * @param matrix
 *     the matrix to raise to a given power
 * @param power
 *     the power
 *
 * @return matrix^power
 * @throw Exception if dimensions leeds to undefined terms
 */
Matrix
pow(const Matrix& matrix, int power) throw (Exception)
{
  Matrix powered_matrix;

  if (power > 0)
    {
      powered_matrix = matrix;
    }
  else if (power < 0)
    {
      powered_matrix = matrix.inverse();
    }
  else
    {
      return Matrix::UNITY(matrix.m_rows);
    }

  for (unsigned int i = 1; i < unsigned(abs(power)); i++)
    {
      powered_matrix *= powered_matrix;
    }

  return powered_matrix;
}

/**
 * @param matrix
 *     a matrix
 *
 * @return the determinant of the matrix
 * @throw Exception if this is no square matrix
 */
double
det(const Matrix& matrix) throw (Exception)
{
  return matrix.det();
}

/**
 * @return true if this is a square NxN matrix
 */
bool
Matrix::isSquare() const
{
  if (m_cols == m_rows)
    {
      return true;
    }
  return false;
}

/**
 * @return true if this matrix is invertible (det!=0)
 */
bool
Matrix::isInvertible() const
{
  if (det() == 0)
    {
      return false;
    }
  return true;
}

/**
 * @return true if this matrix is unitary (U*==1/U)
 */
bool
Matrix::isUnitary() const
{
  try
    {
      if (conjugate() == inverse())
        {
          return true;
        }
      else
        {
          return false;
        }
    }
  catch (Exception e)
    {
      return false;
    }
}

/**
 * @param matrix
 *     another matrix to compare dimensions
 *
 * @return true if number of rows and cols are equal in both matrices
 */
bool
Matrix::matchesDimensions(const Matrix& matrix) const
{
  if (m_rows == matrix.m_rows && m_cols == matrix.m_cols)
    {
      return true;
    }
  return false;
}



Matrix
Matrix::UNITY(unsigned int dimension)
{
	Matrix unity_matrix(dimension, dimension);
	for( unsigned int i = 0; i < dimension; i++ ){
		unity_matrix(i, i) = 1;
    }
    return unity_matrix;
}



Matrix
Matrix::ZERO(unsigned int dimension)
{
	Matrix zero_matrix(dimension, dimension);
	for( unsigned int i = 0; i < dimension; i++ ){
		zero_matrix(i, i) = 0;
    }
	return zero_matrix;
}



/**
 * set all elements to 0
 */
void
Matrix::init()
{
  m_values = new double[m_rows * m_cols];

  for (unsigned int i = 0; i < (m_rows * m_cols); i++)
    {
      m_values[i] = 0;
    }
}



void
Matrix::copy(const Matrix& matrix)
{
	double* old_values = m_values;
 	m_values = new double[matrix.m_rows * matrix.m_cols];
  
	if( NULL!=old_values ){
		delete[] old_values;
		old_values = NULL;
	}

	m_rows = matrix.m_rows;
	m_cols = matrix.m_cols;
 	
	for (unsigned int i = 0; i < (m_rows * m_cols); i++){
		m_values[i] = matrix.m_values[i];
    }
}

