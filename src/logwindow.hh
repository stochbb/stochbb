#ifndef LOGWINDOW_HH
#define LOGWINDOW_HH

#include <QMainWindow>
#include <QListWidget>

#include "lib/logger.hh"

class LogModelHandler: public QObject, public stochbb::LogHandlerObj
{
  Q_OBJECT

public:
  LogModelHandler();

  void handleMessage(const stochbb::LogMessage &msg);

signals:
  void message(QString msg);
};


class LogWindow: public QMainWindow
{
  Q_OBJECT

public:
  LogWindow(QWidget *parent=0);

protected slots:
  void onMessage(QString msg);

protected:
  void closeEvent(QCloseEvent *);

protected:
  QListWidget *_messages;
};

#endif // LOGWINDOW_HH
