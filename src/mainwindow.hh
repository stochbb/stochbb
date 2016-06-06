#ifndef MAINWINDOW_HH
#define MAINWINDOW_HH

#include <QMainWindow>
#include "logwindow.hh"

class NetEditWidget;


class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = 0);

protected slots:
  void onNewNetwork();
  void onSave();
  void onSaveAs();
  void onImageExport();
  void onLoad();
  void onCheck();
  void onRun();
  void onHelp();
  void onAbout();
  void updateTitle();
  void onShowLog(bool show);
  void onQuit();

protected:
  NetEditWidget *_netedit;
  LogWindow *_log;
};

#endif // MAINWINDOW_HH
