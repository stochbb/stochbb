#include "logger.hh"
#include <QFileInfo>
#include <locale>

using namespace stochbb;

/* ********************************************************************************************* *
 * Implementation of LogMessage
 * ********************************************************************************************* */
LogMessage::LogMessage()
  : _filename(), _line(0), _level(LDEBUG), _message(), _timestamp(-1)
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

LogMessageStream &
LogMessageStream::operator<<(const Container &obj) {
  obj->print(*this);
  return *this;
}


/* ********************************************************************************************* *
 * Implementation of LogHandlerObj
 * ********************************************************************************************* */
LogHandlerObj::LogHandlerObj(LogMessage::Level level)
  : stochbb::Object(), _minLevel(level)
{
  // pass...
}

LogHandlerObj::~LogHandlerObj() {
  // pass...
}

void
LogHandlerObj::mark() {
  if (isMarked()) { return; }
  stochbb::Object::mark();
}


/* ********************************************************************************************* *
 * Implementation of LogHandler
 * ********************************************************************************************* */
LogHandler::LogHandler(LogHandlerObj *obj)
  : Container(obj), _loghandler(obj)
{
  // pass...
}

LogHandler::LogHandler(const LogHandler &other)
  : Container(other), _loghandler(other._loghandler)
{
  // pass...
}

LogHandler::~LogHandler() {
  // pass...
}

LogHandler &
LogHandler::operator =(const LogHandler &other) {
  stochbb::Container::operator =(other);
  _loghandler = other._loghandler;
  return *this;
}

void
LogHandler::handleMessage(const LogMessage &msg) {
  _loghandler->handleMessage(msg);
}


/* ********************************************************************************************* *
 * Implementation of IOLogHandlerObj
 * ********************************************************************************************* */
IOLogHandlerObj::IOLogHandlerObj(std::ostream &stream, LogMessage::Level level)
  : LogHandlerObj(level), _stream(stream)
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
IOLogHandlerObj::handleMessage(const LogMessage &msg) {
  if (msg.level() < _minLevel) { return; }
  switch (msg.level()) {
  case LogMessage::LDEBUG: _stream << "DEBUG: "; break;
  case LogMessage::LINFO: _stream << "INFO: "; break;
  case LogMessage::LWARNING: _stream << "WARNING: "; break;
  case LogMessage::LERROR: _stream << "ERROR: "; break;
  }
  std::string basename = msg.filename().substr(msg.filename().find_last_of("/\\") + 1);
  _stream << *std::localtime(&msg.timestamp())
          << ", @"  << basename << ":" << msg.linenumber()
          << ": " << msg.message() << "\n";
  _stream.flush();
}


/* ********************************************************************************************* *
 * Implementation of IOLogHandler
 * ********************************************************************************************* */
IOLogHandler::IOLogHandler(std::ostream &stream, LogMessage::Level level)
  : LogHandler(new IOLogHandlerObj(stream, level))
{
  // pass...
}

IOLogHandler::IOLogHandler(const IOLogHandler &other)
  : LogHandler(other)
{
  // pass...
}

IOLogHandler &
IOLogHandler::operator =(const IOLogHandler &other) {
  LogHandler::operator =(other);
  return *this;
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
  std::list<LogHandlerObj *>::iterator handler = _handler.begin();
  for (; handler != _handler.end(); handler++) {
    (*handler)->unref();
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
  std::list<LogHandlerObj *>::iterator handler = self->_handler.begin();
  for (; handler != self->_handler.end(); handler++) {
    (*handler)->handleMessage(msg);
  }
}

void
Logger::addHandler(const LogHandler &handler) {
  handler->ref();
  Logger().get()->_handler.push_back(
        reinterpret_cast<LogHandlerObj *>(*handler));
}
