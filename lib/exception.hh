#ifndef __SBB_EXCEPTION_HH__
#define __SBB_EXCEPTION_HH__

#include <exception>
#include <sstream>


namespace sbb {

/** Base class of all exceptions. */
class Error: public std::exception, public std::stringstream
{
public:
  /** Constructor. */
  Error();
  /** Copy constructor. */
  Error(const Error &other);
  /** Destructor. */
  virtual ~Error() throw ();
  /** Returns the error message. */
  virtual const char *what() const throw ();
};


/** This exception gets thrown if the parser fails. */
class ParserError: public Error
{
public:
  /** Constructor. */
  ParserError();
  /** Copy constructor. */
  ParserError(const ParserError &other);
  /** Destructor. */
  virtual ~ParserError() throw ();
};

/** This exception gets thrown if an assumption (usually about independence)
 * is not met. */
class AssumptionError: public Error
{
public:
  /** Constructor. */
  AssumptionError();
  /** Copy constructor. */
  AssumptionError(const AssumptionError &other);
  /** Destructor. */
  virtual ~AssumptionError() throw ();
};

}

#endif
