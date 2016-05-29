#ifndef VECTOR_H_
#define VECTOR_H_

#include <iostream>
#include <vector>
#include <cmath>

#include "Exception.hpp"

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

	friend ostream&
	operator<<(ostream& os, const Vector& vector);

	double&
	operator[](unsigned int i) throw (Exception);
	
	const double&
	at(unsigned int i) const throw (Exception);

	Vector&
	operator=(const Vector &vector);

	bool
	operator==(const Vector& vector) const;

	bool
	operator!=(const Vector& vector) const;

	Vector
	operator+(const Vector& vector) const throw (Exception);

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
	operator*(const Vector& vector) const throw (Exception);

	Vector
	operator/(double divisor) const throw (Exception);

	unsigned int
  	size() const;

	unsigned int
  	dim() const;

	double
	norm() const;

	friend double
	abs(const Vector& vector);

	Vector
	normalize() const throw (Exception);

	friend Vector
	normalize(const Vector& vector) throw (Exception);

	double
	angle(const Vector& v) const throw (Exception);

	friend double
	angle(const Vector& a, const Vector& b) throw (Exception);

	Vector
	cross(const Vector& v) const throw (Exception);

	friend Vector
	cross(const Vector& a, const Vector& b) throw (Exception);

	bool
	isNullvector() const;

private:

	vector<double> m_values;
};

#endif /* VECTOR_H_ */
