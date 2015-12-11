#ifndef LOGGER_H
#define LOGGER_H

#include <sstream>
#include <string>
#include <ctime>
#include <list>

#include "object.hh"

/** A log message. */
class LogMessage
{
public:
  /** Specifies possible log levels. */
  typedef enum {
    DEBUG = 0,  ///< Debug information.
    INFO,       ///< Runtime & status information.
    WARNING,    ///< Warnings and minor issues.
    ERROR       ///< Errors and other major issues.
  } Level;

public:
  /** Empty constructor. */
  LogMessage();
  /** Constructor from @c filename, @c line, log @c level and @c message text. */
  LogMessage(const std::string &filename, int line, Level level, const std::string &message);
  /** Copy constructor. */
  LogMessage(const LogMessage &other);
  /** Assignment operator. */
  LogMessage &operator=(const LogMessage &other);

  /** Returns the filename, where the message originated. */
  const std::string &filename() const;
  /** Returns the line number, where the message originated. */
  int linenumber() const;
  /** Returns the level of the message. */
  Level level() const;
  /** Returns the actual message. */
  const std::string &message() const;
  /** Returns the timestamp of the message. */
  const std::time_t &timestamp() const;

protected:
  /** The filename where the message originated. */
  std::string _filename;
  /** The line number where the message originated. */
  int _line;
  /** The message level. */
  Level _level;
  /** The actual message. */
  std::string _message;
  /** The timestamp of the message. */
  std::time_t _timestamp;
};


/** A stream object assembling a log message. Upon destruction, a log message gets assembled and
 * passed to the logger. */
class LogMessageStream: public std::stringstream
{
public:
  /** Constructs a ne log message stream at the position @c filename & @c line with the specified
   * log @c level. Usually the @c logDebug, @c logInfo, @c logWarning or @c logError macros are
   * used to instantiate a @c LogMessageStream. */
  LogMessageStream(const std::string &filename, int line, LogMessage::Level level);
  /** Copy constructor. */
  LogMessageStream(const LogMessageStream &other);
  /** Destructor. Upon destruction, a log message will be assembled and send to the @c Logger
   * instance. */
  virtual ~LogMessageStream();

protected:
  /** The name of the file where the message originated. */
  std::string _filename;
  /** The line where the message originated. */
  int _line;
  /** Level of the log message. */
  LogMessage::Level _level;
};


/** The base class of all log-handlers. */
class LogHandlerObj
{
protected:
  /** Hidden constructor. */
  LogHandlerObj(LogMessage::Level level);

public:
  /** Destructor. */
  virtual ~LogHandlerObj();

  /** Needs to be implemented to handle log messages. */
  virtual void handleMessage(const LogMessage &msg) = 0;

protected:
  /** Minimum log level to process. */
  LogMessage::Level _minLevel;
};


/** Serializes log messages to the given file. */
class IOLogHandlerObj: public LogHandlerObj
{
public:
  /** Constructor.
   * @param level Spicifies the mimimum log level to process.
   * @param stream Specifies the output stream. */
  IOLogHandlerObj(std::ostream &stream, LogMessage::Level level=LogMessage::DEBUG);

  /** Implements the @c LogHandler interface. */
  void handleMessage(const LogMessage &msg);

protected:
  /** A textstream to serialize into. */
  std::ostream &_stream;
};


/** A singleton logger class. */
class Logger
{
protected:
  /** Hidden constructor. */
  Logger();

public:
  /** Destructor. */
  virtual ~Logger();

  /** Logs a message. */
  static void log(const LogMessage &msg);
  /** Adds a handler to the logger, the ownership is transferred to the @c Logger. */
  static void addHandler(LogHandlerObj *handler);

protected:
  /** Factory method. */
  static Logger *get();

protected:
  /** The list of registered handlers. */
  std::list<LogHandlerObj *> _handler;

protected:
  /** The singleton instance. */
  static Logger *_instance;
};

/** Convenience macro to create a @c LogMessageStream with log level "DEBUG". */
#define logDebug()   (LogMessageStream(__FILE__, __LINE__, LogMessage::DEBUG))
/** Convenience macro to create a @c LogMessageStream with log level "INFO". */
#define logInfo()    (LogMessageStream(__FILE__, __LINE__, LogMessage::INFO))
/** Convenience macro to create a @c LogMessageStream with log level "Warning". */
#define logWarning() (LogMessageStream(__FILE__, __LINE__, LogMessage::WARNING))
/** Convenience macro to create a @c LogMessageStream with log level "ERROR". */
#define logError()   (LogMessageStream(__FILE__, __LINE__, LogMessage::ERROR))

#endif // LOGGER_H
