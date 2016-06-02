#ifndef MAINWINDOW_HH
#define MAINWINDOW_HH

#include <QMainWindow>

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
  void onLoad();
  void onCheck();
  void onRun();

protected:
  NetEditWidget *_netedit;
};

#endif // MAINWINDOW_HH
