#include "exception.hh"

using namespace stochbb;


Error::Error()
  : std::exception(), std::stringstream()
{
  this->str("");
}

Error::Error(const Error &other)
  : std::exception(), std::stringstream()
{
  (*this) << other.str();
}

Error::~Error() throw()  {
  // pass...
}

const char *
Error::what() const throw() {
  return this->str().c_str();
}


ParserError::ParserError()
  : Error()
{
  // pass...
}

ParserError::ParserError(const ParserError &other)
  : Error(other)
{
  // pass...
}

ParserError::~ParserError() throw () {
  // pass...
}


AssumptionError::AssumptionError()
  : Error()
{
  // pass...
}

AssumptionError::AssumptionError(const AssumptionError &other)
  : Error(other)
{
  // pass...
}

AssumptionError::~AssumptionError() throw () {
  // pass...
}


TypeError::TypeError()
  : Error()
{
  // pass...
}

TypeError::TypeError(const TypeError &other)
  : Error(other)
{
  // pass...
}

TypeError::~TypeError() throw () {
  // pass...
}
