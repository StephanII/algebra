#include "Matrix.hpp"

#include <sstream>



Matrix::Matrix() :
	m_row_vectors(3)
{ }



Matrix::Matrix(unsigned int rows, unsigned int cols) throw (Exception) :
	m_row_vectors(rows,Vector(cols))
{ 
	if (rows <= 0 || cols <= 0){
    	stringstream ss; ss << "Constructor error, rows=" << rows << " and columns=" << cols << ", but both must be > 0.";
    	throw Exception(__FILE__, __LINE__, ss.str());
	}
}



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



Matrix
Matrix::operator*(const Matrix& matrix) const throw (Exception)
{
   	unsigned int rows = this->rows();
   	unsigned int cols = this->cols();
   	
	if (cols != matrix.rows()) {
		stringstream ss;
		ss << "Error, matrix multiplication of " << rows << "x" << cols << " and " << matrix.rows() << "x" << matrix.cols()	<< " is not defined";
      	throw Exception(__FILE__, __LINE__, ss.str());
    }

	Matrix product(rows, matrix.cols());
	
	for (unsigned int row = 0; row < rows; row++) {
    	for (unsigned int col = 0; col < matrix.cols(); col++) {
        	for (unsigned int i = 0; i < cols; i++) {
            	product[row][col] += at(row, i) * matrix.at(i, col);
            }
        }
    }

	return product;
}



Vector
Matrix::operator*(const Vector& vector) const throw (Exception)
{
   	unsigned int rows = this->rows();
   	unsigned int cols = this->cols();
   	
	if (cols != vector.dim()) {
    	throw Exception(__FILE__, __LINE__, "Error, vector dimension must be the same as the number of matrix columns");
    }

	Vector product(rows);

	for (unsigned int row = 0; row < rows; row++) {
    	for (unsigned int col = 0; col < cols; col++) {
    		product[row] += at(row, col) * vector.at(col);
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
	if (divisor == 0) {
    	throw Exception(__FILE__, __LINE__, "Error, division by zero");
    }
	return (1 / divisor) * (*this);
}



Matrix
Matrix::operator/(const Matrix& matrix) const throw (Exception)
{
	if (!matrix.isInvertible()) {
		throw Exception(__FILE__, __LINE__, "Error, operand is not invertible");
    }
  	return (*this) * matrix.inverse();
}



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



Matrix
Matrix::submatrix(unsigned int remove_row, unsigned int remove_col) const throw (Exception)
{
   	unsigned int rows = this->rows();
   	unsigned int cols = this->cols();	

	if (remove_col >= cols || remove_row >= rows) {
    	throw Exception(__FILE__, __LINE__, "Error, index out of bounds");
    }

	unsigned int submatrix_rows = rows - 1;
	unsigned int submatrix_cols = cols - 1;
  	Matrix submatrix(submatrix_rows,submatrix_cols);

	for (unsigned int row = 0, src_row = 0; row < submatrix_rows; row++, src_row++) {
    	for (unsigned int col = 0, src_col = 0; col < submatrix_cols; col++, src_col++) {
          	if (src_row == remove_row) {
              	src_row++;
            }
          	if (src_col == remove_col) {
            	src_col++;
            }
			submatrix[row][col] = at(src_row, src_col);
        }
    }

	return submatrix;
}



Matrix
Matrix::transposed() const
{
   	unsigned int rows = this->rows();
   	unsigned int cols = this->cols();
   	
	Matrix transposed_matrix(cols, rows);

	for (unsigned int row = 0; row < rows; row++) {
		for (unsigned int col = 0; col < cols; col++) {
    		transposed_matrix[col][row] = at(row, col);
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
	if (!isSquare()) {
    	throw Exception(__FILE__, __LINE__, "Error, there are no cofactors in non square matrices");
    }
	if (col >= cols() || row >= rows()) {
    	throw Exception(__FILE__, __LINE__, "Error, index out of bounds");
    }

	return pow(double(-1), double(row + col)) * submatrix(row, col).det();
}



Matrix
Matrix::adjugate() const throw (Exception)
{
	if (!isSquare()) {
    	throw Exception(__FILE__, __LINE__, "Error, there are no adjugates of non square matrices");
    }

   	unsigned int rows = this->rows();
   	unsigned int cols = this->cols();
   	
  	Matrix adjugate_matrix(rows, cols);

	for (unsigned int row = 0; row < rows; row++) {
		for (unsigned int col = 0; col < cols; col++) {
          	adjugate_matrix[row][col] = cofactor(col, row);
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
	if (!isSquare()) {
    	throw Exception(__FILE__, __LINE__, "Error, there are no inverse on non square matrices");
    }
	if (!isInvertible()) {
      throw Exception(__FILE__, __LINE__, "Error, matrix is not invertible (determinant=0)");
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
	if (!isSquare()) {
    	stringstream ss;
    	ss << "Error, calculating the determinant of a " << rows() << "x" << cols() << " matrix is not defined";
    	throw Exception(__FILE__, __LINE__, ss.str());
    }
    
   	unsigned int rows = this->rows();
   	unsigned int cols = this->cols();    

	if (rows == 1) {
    	return at(0,0);
    } else if (rows == 2) {
    	return at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0);
    } else if (rows == 3) {
    	return at(0, 0) * at(1, 1) * at(2, 2) + 
    		   at(0, 1) * at(1, 2) * at(2, 0) + 
    		   at(0, 2) * at(1, 0) * at(2, 1) - 
    		   at(0, 2) * at(1, 1) * at(2, 0) - 
    		   at(0, 1) * at(1, 0) * at(2, 2) - 
    		   at(0, 0) * at(1, 2) * at(2, 1);
    } else {
      	double determinant = 0;
      	for (unsigned int col = 0; col < cols; col++) {
		    determinant += at(0, col) * cofactor(0,col);
        }
      	return determinant;
    }
}



double
Matrix::trace() const throw (Exception)
{
	if (!isSquare()) {
    	throw Exception(__FILE__, __LINE__, "Error, trace is only defined on square matrices");
    }

   	unsigned int cols = this->cols();

	double trace = 0;
	for (unsigned int i = 0; i < cols; i++) {
      	trace += at(i, i);
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

	if (power > 0) {
    	powered_matrix = matrix;
    } else if (power < 0) {
    	powered_matrix = matrix.inverse();
    } else {
      	return Matrix::UNITY(matrix.rows());
    }

	for (unsigned int i = 1; i < unsigned(abs(power)); i++) {
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
	if (cols() == rows()) {
    	return true;
    }
  	return false;
}



bool
Matrix::isInvertible() const
{
  	if (det() == 0) {
      	return false;
    }
  	return true;
}



bool
Matrix::isUnitary() const
{
	try {
    	if (conjugate() == inverse()) {
          	return true;
        } else {
			return false;
        }
    } catch (Exception e) {
    	return false;
    }
}



Matrix
Matrix::UNITY(unsigned int dimension)
{
	Matrix unity_matrix(dimension, dimension);
	for( unsigned int i = 0; i < dimension; i++ ){
		unity_matrix[i][i] = 1;
    }
    return unity_matrix;
}



Matrix
Matrix::ZERO(unsigned int dimension)
{
	Matrix zero_matrix(dimension, dimension);
	for( unsigned int i = 0; i < dimension; i++ ){
		zero_matrix[i][i] = 0;
    }
	return zero_matrix;
}


