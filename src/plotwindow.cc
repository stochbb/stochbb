#include "plotwindow.hh"
#include <Eigen/Eigen>
#include <QVBoxLayout>


PlotWindow::PlotWindow(double tmin, double tmax, size_t nstep, const QVector<stochbb::Var> &vars, QWidget *parent)
  : QWidget(parent), _tmin(tmin), _tmax(tmax), _nstep(nstep), _vars(vars)
{
  QCustomPlot *plot = new QCustomPlot();

  double t = _tmin, dt = (_tmax-_tmin)/_nstep;
  Eigen::VectorXd F(_nstep), T(_nstep);
  for (size_t i=0; i<_nstep; i++, t+=dt)
    T(i) = t;
  double ymax = 0;
  for (size_t i=0; i<_vars.size(); i++) {
    _vars[i].density().eval(tmin, tmax, F);
    ymax = std::max(ymax, F.maxCoeff());
    QCPGraph *graph = plot->addGraph();
    for (size_t j=0; j<nstep; j++)
      graph->addData(T[j], F[j]);
  }
  plot->xAxis->setRange(tmin, tmax);
  plot->yAxis->setRange(0, ymax);
  plot->replot();

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(plot);
  layout->setMargin(0);
  setLayout(layout);
}
