#include "plot.hh"
#include <QApplication>
#include "qcustomplot.hh"
#include <iostream>
#include <fstream>


int do_plot(int argc, char *argv[], const Eigen::MatrixXd &out, const std::vector<std::string> &names)
{
  QApplication app(argc, argv);
  QCustomPlot *plot = new QCustomPlot();

  QList<QColor> colors;
  colors << QColor(0, 0, 125) << QColor(125, 0, 0) << QColor(0, 125, 0) << QColor(125, 125, 0)
         << QColor(0, 125, 125) << QColor(125, 0, 125) << QColor(205, 79, 18) << QColor(255, 185, 24)
         << QColor(243, 250, 146) << QColor(105, 151, 102) << QColor(69, 47, 96)
         << QColor(224, 26, 53) << QColor(204, 15, 19) << QColor(63, 61, 153) << QColor(153, 61, 113)
         << QColor(61, 153, 86) << QColor(61, 90, 153) << QColor(153, 61, 144) << QColor(61, 121, 153)
         << QColor(132, 61, 153) << QColor(153, 78, 61) << QColor(98, 153, 61) << QColor(61, 151, 153)
         << QColor(101, 61, 153) << QColor(153, 61, 75);

  size_t N = out.rows(), M = out.cols()-1;
  for (size_t j=0; j<M; j++) {
    QVector<double> x, y; x.resize(N), y.resize(N);
    for (size_t i=0; i<N; i++) { x[i] = out(i,0); y[i] = out(i,j+1); }
    // Add a graph
    plot->addGraph();
    plot->graph(j)->setData(x,y);
    plot->graph(j)->setPen(QPen(colors.at(j % colors.size()), 2));
    if (0 == j) { plot->graph(j)->rescaleAxes(); }
    else { plot->graph(j)->rescaleAxes(true); }
    // Add graph to legend
    plot->graph(j)->setName(names[j].c_str());
    plot->graph(j)->addToLegend();
  }
  plot->legend->setVisible(true);

  // go.
  plot->setMinimumSize(700, 500);
  plot->show();
  app.exec();

  return 0;
}


int
output_csv(Eigen::MatrixXd &out, const std::string &filename) {
  std::ofstream file;
  file.open(filename.c_str());
  if (! file.is_open()) {
    std::cerr << "Cannot open file '" << filename << "' for output." << std::endl;
    return -1;
  }
  output_csv(out, file);
  file.flush(); file.close();
  return 0;
}


int
output_csv(Eigen::MatrixXd &out, std::ostream &stream) {
  for (int i=0; i<out.rows(); i++) {
    stream << out(i,0);
    for (int j=1; j<out.cols(); j++) {
      stream << '\t' << out(i,j);
    }
    stream << std::endl;
  }
  return 0;
}
