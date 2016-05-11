#include "Exception.hpp"

#include <iostream>
#include <sstream>



Exception::Exception(const string &file, unsigned int line, const string &what) throw () :
	m_file(file), m_line(line), m_what(what)
{ }



Exception::Exception(Exception& e, const string &file, unsigned int line, const string &what) throw () :
	Exception(file,line,what + string("\n") + e.desc())
{ 

}



Exception::~Exception() throw ()
{ }



const char*
Exception::what() const throw ()
{
  return desc().c_str();
}



const string
Exception::desc() const throw ()
{
  stringstream line;
  line << m_line;
  return ("Exception in file " + m_file + " at line " + line.str() + " : " + m_what);
}



void
Exception::print() const throw ()
{
  cout << what() << endl;
}
