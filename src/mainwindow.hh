#ifndef MAINWINDOW_HH
#define MAINWINDOW_HH

#include <QMainWindow>
#include <QSettings>
#include "logwindow.hh"
#include <QTranslator>


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
  void onLoad(QAction *action);
  void onCheck();
  void onRun();
  void onHelp();
  void onAbout();
  void updateTitle();
  void onShowLog(bool show);
  void onQuit();

protected:
  void addToRecent(const QString &path);
  void populateRecent();

protected:
  NetEditWidget *_netedit;
  LogWindow *_log;
  QMenu *_recent;
  QSettings _settings;
  QTranslator _translator;
};

#endif // MAINWINDOW_HH
