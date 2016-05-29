#include "Vector.hpp"

#include <sstream>



Vector::Vector() :
	m_values(3)
{
	m_values.at(0) = 0;
	m_values.at(1) = 0;
	m_values.at(2) = 0;		
}



Vector::Vector(unsigned int dimension) :
	m_values(dimension)
{
	for( unsigned int i = 0; i < dimension; i++ ){
		m_values.at(i) = 0.;
    }
}



Vector::Vector(double x, double y) :
	m_values(2)
{
	m_values.at(0) = x;
	m_values.at(1) = y;
}



Vector::Vector(double x, double y, double z) :
	m_values(3)
{
	m_values.at(0) = x;
	m_values.at(1) = y;
	m_values.at(2) = z;
}



Vector::Vector(double x1, double x2, double x3, double x4, double x5, double x6) :
	m_values(6)
{
	m_values.at(0) = x1;
	m_values.at(1) = x2;
	m_values.at(2) = x3;
	m_values.at(3) = x4;
	m_values.at(4) = x5;
	m_values.at(5) = x6;
}



Vector::Vector(const Vector& vector) :
	m_values(vector.m_values)
{ }



Vector::~Vector()
{ }



ostream&
operator<<(ostream& os, const Vector& vector)
{
	if (vector.m_values.size() == 0) {
		return os << "Vector(empty)";
    } else {
      	os << "Vector(" << vector.m_values.at(0);

      	for (unsigned int i = 1; i < vector.m_values.size(); i++) {
          	os << ", " << vector.m_values.at(i);
        }

      	os << ")";
      	return os;
    }
}



double&
Vector::operator[](unsigned int i) throw (Exception)
{
	if (i < m_values.size()) {
    	return m_values.at(i);
    } else {
    	stringstream ss; ss << "Error on operator [], index " << i << " out of range [0," << (m_values.size()-1) << "].";
    	throw Exception(__FILE__, __LINE__, ss.str());
    }
}



const double&
Vector::at(unsigned int i) const throw (Exception)
{
	if (i < m_values.size()) {
    	return m_values.at(i);
    } else {
    	stringstream ss; ss << "Error on element access, index " << i << " out of range [0," << (m_values.size()-1) << "].";
    	throw Exception(__FILE__, __LINE__, ss.str());
    }
}



Vector&
Vector::operator=(const Vector &vector)
{
 	m_values = vector.m_values;
	return *this;
}



bool
Vector::operator==(const Vector& vector) const
{
	if (m_values.size() != vector.m_values.size()) {
		return false;
    }

	for (unsigned int i = 0; i < m_values.size(); i++) {
		if (m_values.at(i) != vector.m_values.at(i)) {
        	return false;
        }
    }

	return true;
}



bool
Vector::operator!=(const Vector& vector) const
{
	return not(*this == vector);
}



Vector
Vector::operator+(const Vector& vector) const throw (Exception)
{
	if (m_values.size() != vector.m_values.size()) {
    	stringstream ss; ss << "Error on operator +, dimension " << m_values.size() << " does not fit " << vector.m_values.size() << ".";
    	throw Exception(__FILE__, __LINE__, ss.str());
    }

	Vector sum(*this);
	
	for (unsigned int i = 0; i < m_values.size(); i++) {
		sum[i] += vector.m_values.at(i);
    }

	return sum;
}



Vector&
Vector::operator+=(const Vector& vector)
{
	*this = *this + vector;
	return *this;
}



Vector
Vector::operator-(const Vector& vector) const
{
	return (*this + (-1 * vector));
}



Vector
operator-(const Vector& vector)
{
	return (-1 * vector);
}



Vector&
Vector::operator-=(const Vector& vector)
{
	*this = *this - vector;
	return *this;
}



Vector
Vector::operator*(double factor) const
{
	Vector product(m_values.size());

	for (unsigned int i = 0; i < m_values.size(); i++) {
    	product[i] = m_values.at(i) * factor;
    }

	return product;
}



Vector
operator*(double factor, const Vector& vector)
{
	Vector product(vector.m_values.size());

	for (unsigned int i = 0; i < vector.m_values.size(); i++) {
    	product[i] = vector.m_values.at(i) * factor;
    }

	return product;
}



double
Vector::operator*(const Vector& vector) const throw (Exception)
{
	if (m_values.size() != vector.m_values.size()) {
    	stringstream ss; ss << "Error on operator *, dimension " << m_values.size() << " does not fit " << vector.m_values.size() << ".";
    	throw Exception(__FILE__, __LINE__, ss.str());
    }
    
	double product = 0.;

  	for (unsigned int i = 0; i < m_values.size(); i++) {
    	product += m_values.at(i) * vector.m_values.at(i);
    }

	return product;
}



Vector
Vector::operator/(double divisor) const throw (Exception)
{
	if (divisor == 0) {
		throw Exception(__FILE__, __LINE__, "Error, division by zero.");
    }

	return (1. / divisor) * (*this);
}



unsigned int
Vector::size() const
{
	return m_values.size();
}


unsigned int
Vector::dim() const
{
	return size();
}



double
Vector::norm() const
{
	double norm = 0;

	for (unsigned int i = 0; i < m_values.size(); i++) {
		norm += m_values[i] * m_values[i];
    }

	return sqrt(norm);
}



double
abs(const Vector& vector)
{
	return vector.norm();
}



Vector
Vector::normalize() const throw (Exception)
{
	if (isNullvector()) {
	    throw Exception(__FILE__, __LINE__, "Error on normalization, nullvector has no direction and can not be normed");
    }
  	return ((*this) / norm());
}



Vector
normalize(const Vector& vector) throw (Exception)
{
	return vector.normalize();
}



double
Vector::angle(const Vector& vector) const throw (Exception)
{
	if (norm() * vector.norm() == 0) {
		throw Exception(__FILE__, __LINE__, "Error, nullvector has no direction ergo no angle.");
    } else {
		double cos_angle = (((*this) * vector) / (this->norm() * vector.norm()));
		if (cos_angle >= 1){ 
			return 0; 
		} else if (cos_angle <= -1){ 
			return M_PI;
		} else {		
	      	return acos(cos_angle);
	    }
    }
}



double
angle(const Vector& a, const Vector& b) throw (Exception)
{
	return a.angle(b);
}



Vector
Vector::cross(const Vector& vector) const throw (Exception)
{
	if (m_values.size() != 3 || vector.m_values.size() != 3) {
      throw Exception(__FILE__, __LINE__, "Error, cross product is only defined in 3 dimensions.");
    } else {
    	Vector cross_product(3);
    	cross_product[0] = m_values[1] * vector.m_values[2] - m_values[2] * vector.m_values[1];
    	cross_product[1] = m_values[2] * vector.m_values[0] - m_values[0] * vector.m_values[2];
      	cross_product[2] = m_values[0] * vector.m_values[1] - m_values[1] * vector.m_values[0];
    	return cross_product;
    }
}



Vector
cross(const Vector& a, const Vector& b) throw (Exception)
{
	return a.cross(b);
}



bool
Vector::isNullvector() const
{
	for (unsigned int i = 0; i < dim(); i++) {
    	if (m_values[i] != 0) {
    		return false;
        }
    }
	return true;
}


