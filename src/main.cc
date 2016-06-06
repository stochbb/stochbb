#include "lib/api.hh"
#include "lib/option_parser.hh"
#include "lib/logger.hh"
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
