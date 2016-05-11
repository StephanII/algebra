/*
 * Vector.h
 *
 *  Created on: 02.09.2009
 *      Author: stephanreimann
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <iostream>
#include <cmath>

#include "Exception.hpp"
#include "Matrix.hpp"

using namespace std;

class Vector
{
public:

	Vector();

	Vector(unsigned int dimension);

	Vector(double x, double y);

	Vector(double x, double y, double z);
  
	Vector(double x1, double x2, double x3, double x4, double x5, double x6); 

	Vector(const Vector& vector);

	virtual
	~Vector();

	double&
	operator()(unsigned int i) throw (Exception);

	double&
	operator[](unsigned int i) throw (Exception);

	Vector&
	operator=(const Vector &vector);

	bool
	operator==(const Vector& vector) const;

	bool
	operator!=(const Vector& vector) const;

	Vector
	operator+(const Vector& vector) const;

  	Vector&
	operator+=(const Vector& vector);

	Vector
	operator-(const Vector& vector) const;

	friend Vector
	operator-(const Vector& vector);

	Vector&
	operator-=(const Vector& vector);

	Vector
	operator*(double factor) const;

	friend Vector
	operator*(double factor, const Vector& vector);

  	double
	operator*(const Vector& vector) const;

	Vector
	operator/(double divisor) const throw (Exception);

	friend ostream&
	operator<<(ostream& os, const Vector& vector);

	unsigned int
  	getDimension() const;

	double
	get(unsigned int i) const throw (Exception);

	void
	set(unsigned int i, double value) throw (Exception);

	double
	norm() const;

	friend double
	abs(const Vector& vector);

	Vector
	normalized() const throw (Exception);

	friend Vector
	normalized(const Vector& vector) throw (Exception);

	double
	angle(const Vector& v) const throw (Exception);

	friend double
	angle(const Vector& a, const Vector& b) throw (Exception);

	Vector
	crossProduct(const Vector& v) const throw (Exception);

	friend Vector
	crossProduct(const Vector& a, const Vector& b) throw (Exception);

	bool
	isNullvector() const;

private:

	void
	copy(const Vector& vector);

	unsigned int m_dimension;

	double* m_values;
};

#endif /* VECTOR_H_ */
