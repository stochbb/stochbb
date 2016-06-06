#include <stochbb/api.hh>
#include <stochbb/option_parser.hh>
#include <stochbb/logger.hh>
#include "qcustomplot.hh"
#include "mainwindow.hh"

#include <iostream>
#include <fstream>

#include <QApplication>


using namespace stochbb;

int main(int argc, char *argv[]) {

  //stochbb::Logger::addHandler(stochbb::IOLogHandler());

  QApplication app(argc, argv);

  MainWindow *win = new MainWindow();
  win->show();

  app.exec();

  return 0;
}
