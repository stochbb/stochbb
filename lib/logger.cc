#include "logger.hh"
#include <QFileInfo>
#include <locale>

/* ********************************************************************************************* *
 * Implementation of LogMessage
 * ********************************************************************************************* */
LogMessage::LogMessage()
  : _filename(), _line(0), _level(DEBUG), _message(), _timestamp(-1)
{
  // pass...
}

LogMessage::LogMessage(const std::string &filename, int line, Level level, const std::string &message)
  : _filename(filename), _line(line), _level(level), _message(message),
    _timestamp(std::time(0))
{
  // pass...
}

LogMessage::LogMessage(const LogMessage &other)
  : _filename(other._filename), _line(other._line), _level(other._level), _message(other._message),
    _timestamp(other._timestamp)
{
  // pass...
}

LogMessage &
LogMessage::operator =(const LogMessage &other) {
  _filename  = other._filename;
  _line      = other._line;
  _level     = other._level;
  _timestamp = other._timestamp;
  return *this;
}

const std::string &
LogMessage::filename() const {
  return _filename;
}

int
LogMessage::linenumber() const {
  return _line;
}

LogMessage::Level
LogMessage::level() const {
  return _level;
}

const std::string &
LogMessage::message() const {
  return _message;
}

const std::time_t &
LogMessage::timestamp() const {
  return _timestamp;
}


/* ********************************************************************************************* *
 * Implementation of LogMessageStream
 * ********************************************************************************************* */
LogMessageStream::LogMessageStream(const std::string &filename, int line, LogMessage::Level level)
  : std::stringstream(), _filename(filename), _line(line), _level(level)
{
  // pass...
}

LogMessageStream::LogMessageStream(const LogMessageStream &other)
  : std::stringstream(), _filename(other._filename), _line(other._line), _level(other._level)
{
  this->str(other.str());
}

LogMessageStream::~LogMessageStream() {
  Logger::log(LogMessage(_filename, _line, _level, this->str()));
}


/* ********************************************************************************************* *
 * Implementation of LogHandler
 * ********************************************************************************************* */
LogHandler::LogHandler(LogMessage::Level level)
  : _minLevel(level)
{
  // pass...
}

LogHandler::~LogHandler() {
  // pass...
}


/* ********************************************************************************************* *
 * Implementation of IOLogHandler
 * ********************************************************************************************* */
IOLogHandler::IOLogHandler(std::ostream &stream, LogMessage::Level level)
  : LogHandler(level), _stream(stream)
{
  // pass...
}

std::ostream & operator<<(std::ostream &stream, const std::tm &time) {
  char datetime[257];
  size_t len = std::strftime(datetime, 256, "%c", &time);
  datetime[len] = 0;
  stream << datetime;
  return stream;
}

void
IOLogHandler::handleMessage(const LogMessage &msg) {
  if (msg.level() < _minLevel) { return; }
  switch (msg.level()) {
  case LogMessage::DEBUG: _stream << "DEBUG: "; break;
  case LogMessage::INFO: _stream << "INFO: "; break;
  case LogMessage::WARNING: _stream << "WARNING: "; break;
  case LogMessage::ERROR: _stream << "ERROR: "; break;
  }
  std::string basename = msg.filename().substr(msg.filename().find_last_of("/\\") + 1);
  _stream << *std::localtime(&msg.timestamp())
          << ", @"  << basename << ":" << msg.linenumber()
          << ": " << msg.message() << "\n";
  _stream.flush();
}


/* ********************************************************************************************* *
 * Implementation of Logger
 * ********************************************************************************************* */
Logger *Logger::_instance = 0;

Logger::Logger()
  : _handler()
{
  // pass...
}

Logger::~Logger() {
  std::list<LogHandler *>::iterator handler = _handler.begin();
  for (; handler != _handler.end(); handler++) {
    delete (*handler);
  }
}

Logger *
Logger::get() {
  if (0 == _instance) {
    _instance = new Logger();
  }
  return _instance;
}

void
Logger::log(const LogMessage &msg) {
  Logger *self = Logger().get();
  std::list<LogHandler *>::iterator handler = self->_handler.begin();
  for (; handler != self->_handler.end(); handler++) {
    (*handler)->handleMessage(msg);
  }
}

void
Logger::addHandler(LogHandler *handler) {
  Logger *self = Logger().get();
  self->_handler.push_back(handler);
}
