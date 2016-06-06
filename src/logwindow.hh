#ifndef LOGWINDOW_HH
#define LOGWINDOW_HH

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


class LogWindow: public QListWidget
{
  Q_OBJECT

public:
  LogWindow(QWidget *parent=0);

protected slots:
  void onMessage(QString msg);
};

#endif // LOGWINDOW_HH
