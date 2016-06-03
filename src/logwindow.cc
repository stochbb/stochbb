#include "logwindow.hh"
#include <locale>


LogModelHandler::LogModelHandler()
  : QObject(), LogHandlerObj(stochbb::LogMessage::LDEBUG)
{
  // pass...
}

void
LogModelHandler::handleMessage(const stochbb::LogMessage &msg) {
  std::stringstream stream;
  switch (msg.level()) {
  case stochbb::LogMessage::LDEBUG: stream << "DEBUG: "; break;
  case stochbb::LogMessage::LINFO: stream << "INFO: "; break;
  case stochbb::LogMessage::LWARNING: stream << "WARNING: "; break;
  case stochbb::LogMessage::LERROR: stream << "ERROR: "; break;
  }
  std::string basename = msg.filename().substr(msg.filename().find_last_of("/\\") + 1);
  stream << ", @"  << basename << ":" << msg.linenumber()
         << ": " << msg.message() << "\n";
  emit message(QString::fromStdString(stream.str()));
}


LogWindow::LogWindow(QWidget *parent)
  : QMainWindow(parent), _messages(0)
{
  _messages = new QListWidget();

  LogModelHandler *handler = new LogModelHandler();
  connect(handler, SIGNAL(message(QString)), this, SLOT(onMessage(QString)));
  stochbb::Logger::addHandler(handler);

  setCentralWidget(_messages);
}

void
LogWindow::onMessage(QString msg) {
  _messages->insertItem(0, msg);
}

void
LogWindow::closeEvent(QCloseEvent *evt) {
  this->hide();
}
