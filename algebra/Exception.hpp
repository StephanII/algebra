/*
 * Exception.h
 *
 *  Created on: 27.08.2009
 *      Author: stephanreimann
 */

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include <exception>
#include <string>
#include <sstream>

using namespace std;

template <typename T>
  string NTS ( T number )
  {
     ostringstream ss;
     ss << number;
     return ss.str();
  }

class Exception : public exception
{

public:
  Exception(const string &file, unsigned int line, const string &what) throw ();
  
  Exception(Exception& e, const string &file, unsigned int line, const string &what) throw ();

  virtual
  ~Exception() throw ();

  const char*
  what() const throw ();

  const string
  desc() const throw ();

  void
  print() const throw ();

private:

  string m_file;

  unsigned int m_line;

  string m_what;
};

#endif /* EXCEPTION_H_ */
