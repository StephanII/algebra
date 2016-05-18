#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <cmath>

#include "Exception.hpp"
#include "Vector.hpp"

class Vector;

class Matrix
{
public:

	Matrix();

	Matrix(unsigned int rows, unsigned int cols);

	Matrix(const Matrix& matrix);

	virtual
	~Matrix();

	friend ostream&
	operator<<(ostream& os, const Matrix& matrix);

  	Vector&
  	operator[](unsigned int row) throw (Exception);
  
	const double&
  	at(unsigned int row, unsigned int col) const throw (Exception); 

	Matrix&
	operator=(const Matrix& matrix);

	bool
	operator==(const Matrix& matrix) const;

	bool
	operator!=(const Matrix& matrix) const;

	Matrix
	operator+(const Matrix& matrix) const throw (Exception);

	Matrix&
	operator+=(const Matrix& matrix);

	Matrix
	operator-(const Matrix& matrix) const throw (Exception);

	Matrix&
	operator-=(const Matrix& matrix);

	friend Matrix
	operator-(const Matrix& matrix);

	Matrix
	operator*(double factor) const;

	friend Matrix
	operator*(double factor, const Matrix& matrix);
/*
  Matrix
  operator*(const Matrix& matrix) const throw (Exception);

  Vector
  operator*(const Vector& vector) const throw (Exception);

  Matrix&
  operator*=(const Matrix& matrix) throw (Exception);

  Matrix
  operator/(double divisor) const throw (Exception);

  Matrix
  operator/(const Matrix& matrix) const throw (Exception);

*/
  unsigned int
  rows() const;

  unsigned int
  cols() const;
/*

  Matrix
  submatrix(unsigned int remove_row, unsigned int remove_col) const
      throw (Exception);

  Matrix
  transposed() const;

  friend Matrix
  transposed(const Matrix& matrix);

  double
  cofactor(unsigned int row, unsigned int col) const throw (Exception);

  Matrix
  adjugate() const throw (Exception);

  friend Matrix
  adjugate(const Matrix& matrix) throw (Exception);

  Matrix
  conjugate() const;

  friend Matrix
  conjugate(const Matrix& matrix);

  Matrix
  inverse() const throw (Exception);

  friend Matrix
  inverse(const Matrix& matrix) throw (Exception);

  double
  det() const throw (Exception);

  friend double
  det(const Matrix& matrix) throw (Exception);

  double
  trace() const throw (Exception);

  friend double
  trace(const Matrix& matrix) throw (Exception);

  friend Matrix
  pow(const Matrix& matrix, int power) throw (Exception);

  bool
  isSquare() const;

  bool
  isInvertible() const;

  bool
  isUnitary() const;

  static Matrix
  UNITY(unsigned int dimension);

  static Matrix
  ZERO(unsigned int dimension);
*/
private:

  vector<Vector> m_row_vectors;
};

#endif /* MATRIX_H_ */
