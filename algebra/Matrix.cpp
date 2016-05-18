#include "Matrix.hpp"

#include <sstream>

Matrix::Matrix() :
	m_row_vectors(3)
{ }



Matrix::Matrix(unsigned int rows, unsigned int cols) :
	m_row_vectors(rows,Vector(cols))
{ }



Matrix::Matrix(const Matrix& matrix) :
	m_row_vectors(matrix.m_row_vectors)
{ }



Matrix::~Matrix()
{ }



ostream&
operator<<(ostream& os, const Matrix& matrix)
{
	if (matrix.rows() == 0 || matrix.cols() == 0) {
    	return os << "Matrix(empty)";
    } else {
    	os << "Matrix" << matrix.rows() << "X" << matrix.cols() << "{";
    	for (unsigned int row = 0; row < matrix.rows(); row++) {
        	os << "(" << matrix.at(row,0);
        	for (unsigned int col = 1; col < matrix.cols(); col++) {
	            os << ", " << matrix.at(row,col);
            }
    		os << ")";
        }
   		os << "}";
    	return os;
    }
}



Vector&
Matrix::operator[](unsigned int row) throw (Exception)
{
	if (row >= m_row_vectors.size()) {
    	stringstream ss; ss << "Error on operator [], row " << row << " out of range [0," << m_row_vectors.size() << "].";
    	throw Exception(__FILE__, __LINE__, ss.str());
	} else {
    	return m_row_vectors[row];
    }
}



const double&
Matrix::at(unsigned int row, unsigned int col) const throw (Exception)
{
	if (row >= m_row_vectors.size()) {
    	stringstream ss; ss << "Error on Matrix element access, index " << row << "," << col << " out of range " << this->rows() << "x" << this->cols() << ".";
    	throw Exception(__FILE__, __LINE__, ss.str());
	} else {
    	return m_row_vectors.at(row).at(col);
    }
}



Matrix&
Matrix::operator=(const Matrix &matrix)
{
	m_row_vectors = matrix.m_row_vectors;
	return *this;
}



bool
Matrix::operator==(const Matrix& matrix) const
{
	if ( not(rows() == matrix.rows() && cols() == matrix.cols()) ) {
		return false;
    }

	for (unsigned int row = 0; row<rows(); row++) {
		for (unsigned int col = 0; col<cols(); col++) {	
			if (at(row,col) != matrix.at(row,col)) {
	        	return false;
	        }
        }
    }

	return true;
}



bool
Matrix::operator!=(const Matrix& matrix) const
{
  return not (*this == matrix);
}



Matrix
Matrix::operator+(const Matrix& matrix) const throw (Exception)
{
   	unsigned int rows = this->rows();
   	unsigned int cols = this->cols();

	if ( not(rows == matrix.rows() && cols == matrix.cols()) ) {
    	stringstream ss; ss << "Error on Matrix operator+, dimension " << rows << "x" << cols << " does not match " << matrix.rows() << "x" << matrix.cols() << ".";
    	throw Exception(__FILE__, __LINE__, ss.str());
    } else {
		Matrix sum(rows, cols);
	    for (unsigned int row = 0; row < rows; row++) {
			for (unsigned int col = 0; col < cols; col++) {
        		sum[row][col] = at(row,col) + matrix.at(row,col);
        	}
    	}
    	return sum;
    }
}



Matrix&
Matrix::operator+=(const Matrix& matrix) 
{
	(*this) = (*this) + matrix;
	return (*this);
}



Matrix
Matrix::operator-(const Matrix& matrix) const throw (Exception)
{
	return (*this) + (-1 * matrix);
}



Matrix&
Matrix::operator-=(const Matrix& matrix)
{
	(*this) = (*this) - matrix;
	return (*this);
}



Matrix
operator-(const Matrix& matrix)
{
	return (-1 * matrix);
}



Matrix
Matrix::operator*(double factor) const
{
   	unsigned int rows = this->rows();
   	unsigned int cols = this->cols();
   	
	Matrix product(rows, cols);

	for (unsigned int row = 0; row < rows; row++) {
    	for (unsigned int col = 0; col < cols; col++) {
        	product[row][col] = factor * at(row,col);
        }
    }

	return product;
}



Matrix
operator*(double factor, const Matrix& matrix)
{
	return matrix * factor;
}


/*
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



Matrix&
Matrix::operator*=(const Matrix& matrix) throw (Exception)
{
  (*this) = (*this) * matrix;
  return (*this);
}



Matrix
Matrix::operator/(double divisor) const throw (Exception)
{
  if (divisor == 0)
    {
      throw Exception(__FILE__, __LINE__, "Error, division by zero");
    }
  return (1 / divisor) * (*this);
}



Matrix
Matrix::operator/(const Matrix& matrix) const throw (Exception)
{
  if (!matrix.isInvertible())
    {
      throw Exception(__FILE__, __LINE__, "Error, operand is not invertible");
    }
  return (*this) * matrix.inverse();
}



*/
unsigned int
Matrix::rows() const
{
	return m_row_vectors.size();
}



unsigned int
Matrix::cols() const
{
	return m_row_vectors.at(0).size();
}

/*



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



Matrix
transposed(const Matrix& matrix)
{
  return matrix.transposed();
}



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



Matrix
adjugate(const Matrix& matrix) throw (Exception)
{
  return matrix.adjugate();
}



Matrix
Matrix::conjugate() const
{
  return transposed();
}



Matrix
conjugate(const Matrix& matrix)
{
  return matrix.conjugate();
}



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



Matrix
inverse(const Matrix& matrix) throw (Exception)
{
  return matrix.inverse();
}



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



double
trace(const Matrix& matrix) throw (Exception)
{
  return matrix.trace();
}



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



double
det(const Matrix& matrix) throw (Exception)
{
  return matrix.det();
}



bool
Matrix::isSquare() const
{
  if (m_cols == m_rows)
    {
      return true;
    }
  return false;
}



bool
Matrix::isInvertible() const
{
  if (det() == 0)
    {
      return false;
    }
  return true;
}



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

*/
