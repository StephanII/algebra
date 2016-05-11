#include "Vector.hpp"



Vector::Vector() :
	m_dimension(3)
{
	m_values = new double[3];
	m_values[0] = 0;
	m_values[1] = 0;
	m_values[2] = 0;
}

/**
 * standard constructor
 * 
 * @param dimension 
 *     the vectors dimension
 */
Vector::Vector(unsigned int dimension) :
  m_dimension(dimension)
{
	m_values = new double[dimension];
	for( unsigned int i = 0; i < dimension; i++ ){
		m_values[i] = 0.;
    }
}

/**
 * constructor for 2 dimensions
 *
 * @param x 
 *     the x coordinate
 * @param y 
 *     the y coordinate
 */
Vector::Vector(double x, double y) :
  m_dimension(2)
{
  m_values = new double[2];
  m_values[0] = x;
  m_values[1] = y;
}

/**
 * constructor for 3 dimensions
 *
 * @param x 
 *     the x coordinate
 * @param y 
 *     the y coordinate
 * @param z 
 *     the z coordinate
 */
Vector::Vector(double x, double y, double z) :
  m_dimension(3)
{
  m_values = new double[3];
  m_values[0] = x;
  m_values[1] = y;
  m_values[2] = z;
}



Vector::Vector(double x1, double x2, double x3, double x4, double x5, double x6) :
	m_dimension(6)
{
	m_values = new double[6];
	m_values[0] = x1;
	m_values[1] = x2;
	m_values[2] = x3;
	m_values[3] = x4;
	m_values[4] = x5;
	m_values[5] = x6;
}



Vector::Vector(const Vector& vector) :
	m_dimension(vector.m_dimension),
	m_values(NULL)
{
	copy(vector);
}



/**
 * destructor
 */
Vector::~Vector()
{
  delete[] m_values;
}

/**
 * @param i 
 *     the index
 * 
 * @return the i-th element of this vector 
 * 
 * @throw Exception
 *     if index is out of range
 */
double&
Vector::operator()(unsigned int i) throw (Exception)
{
  if (i < m_dimension)
    {
      return *(m_values + i);
    }
  else
    {
      throw Exception(__FILE__, __LINE__, "Error, index out of range");
    }
}

/**
 * @param i
 *     the index
 *
 * @return the i-th element of this vector
 *
 * @throw Exception
 *     if index is out of range
 */
double&
Vector::operator[](unsigned int i) throw (Exception)
{
  if (i < m_dimension)
    {
      return *(m_values + i);
    }
  else
    {
      throw Exception(__FILE__, __LINE__, "Error, index out of range");
    }
}



Vector&
Vector::operator=(const Vector &vector)
{
	copy(vector);
	return *this;
}



/**
 * @param vector
 *     a vector to compare
 * 
 * @return true, if bot vectors are equal
 */
bool
Vector::operator==(const Vector& vector) const
{
  if (m_dimension != vector.m_dimension)
    {
      return false;
    }

  for (unsigned int i = 0; i < m_dimension; i++)
    {
      if (m_values[i] != vector.m_values[i])
        {
          return false;
        }
    }

  return true;
}

/**
 * @param vector
 *     a vector to compare
 * 
 * @return false, if bot vectors are equal
 */
bool
Vector::operator!=(const Vector& vector) const
{
  return !(*this == vector);
}

/**
 * @param vector
 *     a vector to add
 *
 * @return the sum of this and the argument
 */
Vector
Vector::operator+(const Vector& vector) const
{
  unsigned int min_dimension = m_dimension <= vector.m_dimension ? m_dimension
      : vector.m_dimension;

  Vector sum(m_dimension >= vector.m_dimension ? *this : vector);

  for (unsigned int i = 0; i < min_dimension; i++)
    {
      sum[i] = m_values[i] + vector.m_values[i];
    }

  return sum;
}

/**
 * @param vector
 *     a vector to add
 * 
 * @return *this after the vector was added
 */
Vector&
Vector::operator+=(const Vector& vector)
{
  *this = *this + vector;
  return *this;
}

/**
 * @param vector
 *     a vector to substract
 *
 * @return the difference of this and the argument
 */
Vector
Vector::operator-(const Vector& vector) const
{
  return (*this + (-1 * vector));
}

/**
 * @return a vector where all components have inverse signs
 */
Vector
operator-(const Vector& vector)
{
  return (-1 * vector);
}

/**
 * @param vector
 *     a vector to substract
 * 
 * @return *this after the vector was substracted
 */
Vector&
Vector::operator-=(const Vector& vector)
{
  *this = *this - vector;
  return *this;
}

/**
 * @param factor
 *     a factor to multiply
 *
 * @return a vector stretched by the given factor
 */
Vector
Vector::operator*(double factor) const
{
  Vector product(m_dimension);

  for (unsigned int i = 0; i < m_dimension; i++)
    {
      product[i] = m_values[i] * factor;
    }

  return product;
}

/**
 * @param factor
 *     a factor to multiply
 * @param vector
 *     a vector to multiply
 *
 * @return the vector stretched by the given factor
 */
Vector
operator*(double factor, const Vector& vector)
{
  Vector product(vector.m_dimension);

  for (unsigned int i = 0; i < vector.m_dimension; i++)
    {
      product[i] = vector.m_values[i] * factor;
    }

  return product;
}

/**
 * the scalar product
 * 
 * @param vector
 *     a vector to multiply
 *
 * @return the scalar product x1*x2 + y1*y2 + z1*z2 + ...
 */
double
Vector::operator*(const Vector& vector) const
{
  double product = 0.;

  unsigned int min_dimension = m_dimension <= vector.m_dimension ? m_dimension
      : vector.m_dimension;

  for (unsigned int i = 0; i < min_dimension; i++)
    {
      product += m_values[i] * vector.m_values[i];
    }

  return product;
}

/**
 * @param divisor
 *     the divisor
 *
 * @return the matrix where all elements are divided by the divisor
 * @throw Exception if the divisor = 0 (division by zero)
 */
Vector
Vector::operator/(double divisor) const throw (Exception)
{
  if (divisor == 0)
    {
      throw Exception(__FILE__, __LINE__, "Error, division by zero");
    }
  return (1 / divisor) * (*this);
}

/**
 * @param os 
 *     an ofstream object
 * @param vector
 *     a vector to stream
 *
 * @return the ofstream object containing the vector data
 */
ostream&
operator<<(ostream& os, const Vector& vector)
{
  if (vector.m_dimension == 0)
    {
      return os << "Vector0D(empty)";
    }
  else if (vector.m_dimension == 1)
    {
      return os << "Vector1D(" << vector.m_values[0] << ")";
    }
  else
    {
      os << "Vector" << vector.m_dimension << "D(" << vector.m_values[0];

      for (unsigned int i = 1; i < vector.m_dimension; i++)
        {
          os << ", " << vector.m_values[i];
        }

      os << ")";
      return os;
    }
}

/**
 * @return the number of coordinates, the vector represents
 */
unsigned int
Vector::getDimension() const
{
  return m_dimension;
}

/**
 * @param i
 *     the index
 *
 * @return the value at index i
 * @throw Exception if index is out of bounds
 */
double
Vector::get(unsigned int i) const throw (Exception)
{
  if (i < m_dimension)
    {
      return m_values[i];
    }
  else
    {
      throw Exception(__FILE__, __LINE__, "Error, index out of range");
    }
}

/**
 * @param i
 *     the index
 * @param value
 *     the value to set
 *
 * @throw Exception if index is out of bounds
 */
void
Vector::set(unsigned int i, double value) throw (Exception)
{
  if (i < m_dimension)
    {
      m_values[i] = value;
    }
  else
    {
      throw Exception(__FILE__, __LINE__, "Error, index out of range");
    }
}

/**
 * @return the length of this vector
 */
double
Vector::norm() const
{
  double norm = 0;

  for (unsigned int i = 0; i < m_dimension; i++)
    {
      norm += m_values[i] * m_values[i];
    }

  return sqrt(norm);
}

/**
 * @param vector
 *     a vector
 *
 * @return the length of the vector
 */
double
abs(const Vector& vector)
{
  return vector.norm();
}

/**
 * @return the normed vector (same direction but length=1)
 * @throw Exception if this is a nullvector
 */
Vector
Vector::normalized() const throw (Exception)
{
  if (isNullvector())
    {
      throw Exception(__FILE__, __LINE__,
          "Error, nullvector has no direction and can not be normed");
    }
  return ((*this) / norm());
}

/**
 * @param vector
 *     a vector to normalize
 *
 * @return the normed vector (same direction but length=1)
 * @throw Exception if this is a nullvector
 */
Vector
normalized(const Vector& vector) throw (Exception)
{
  return vector.normalized();
}

/**
 * @param v
 *     another vector
 *
 * @return the included angle between both vectors this*v/|this|*|v|
 */
double
Vector::angle(const Vector& v) const throw (Exception)
{
  if (norm() * v.norm() == 0)
    {
      throw Exception(__FILE__, __LINE__,
          "Error, nullvector has no direction ergo no angle");
    }
  else
    {
      return ((*this * v) / this->norm() * v.norm());
    }
}

/**
 * @param a
 *     first vector
 * @param b
 *     second vector
 *
 * @return the included angle between both vectors a*b/|a|*|b|
 */
double
angle(const Vector& a, const Vector& b) throw (Exception)
{
  return a.angle(b);
}

/**
 * @param v
 *     another vector
 *
 * @return the cross product of both vectors
 * @throw Exception if one ore both vectors have not 3 dimensions
 */
Vector
Vector::crossProduct(const Vector& v) const throw (Exception)
{
  if (m_dimension != 3 || v.m_dimension != 3)
    {
      throw Exception(__FILE__, __LINE__,
          "Error, cross product is only defined on 3 dimensional space");
    }
  else
    {
      Vector cross_product(3);
      cross_product[0] = m_values[2] * v.m_values[3] - m_values[3]
          * v.m_values[2];
      cross_product[1] = m_values[3] * v.m_values[1] - m_values[1]
          * v.m_values[3];
      cross_product[2] = m_values[1] * v.m_values[2] - m_values[2]
          * v.m_values[1];
      return cross_product;
    }
}

/**
 * @param a
 *     first vector
 * @param b
 *     second vector
 *
 * @return the cross product of both vectors
 * @throw Exception if one ore both vectors have not 3 dimensions
 */
Vector
crossProduct(const Vector& a, const Vector& b) throw (Exception)
{
  return a.crossProduct(b);
}

/**
 * @return true if this is a nullvector
 */
bool
Vector::isNullvector() const
{
  for (unsigned int i = 0; i < m_dimension; i++)
    {
      if (m_values[i] != 0)
        {
          return false;
        }
    }
  return true;
}



void
Vector::copy(const Vector& vector)
{
	double* old_values = m_values;
 	m_values = new double[vector.m_dimension];
  
	if( NULL!=old_values ){
		delete[] old_values;
		old_values = NULL;
	}

	m_dimension = vector.m_dimension;
 	
	for (unsigned int i = 0; i < m_dimension; i++){
		m_values[i] = vector.m_values[i];
    }
}
