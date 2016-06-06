#include "plotwindow.hh"
#include <Eigen/Eigen>

QVector<QColor> colors(
{ QColor(0, 0, 125), QColor(125, 0, 0), QColor(0, 125, 0), QColor(125, 125, 0), QColor(0, 125, 125),
  QColor(125, 0, 125), QColor(205, 79, 18), QColor(255, 185, 24), QColor(243, 250, 146),
  QColor(105, 151, 102), QColor(69, 47, 96), QColor(224, 26, 53), QColor(204, 15, 19),
  QColor(63, 61, 153), QColor(153, 61, 113), QColor(61, 153, 86), QColor(61, 90, 153),
  QColor(153, 61, 144), QColor(61, 121, 153), QColor(132, 61, 153), QColor(153, 78, 61),
  QColor(98, 153, 61), QColor(61, 151, 153), QColor(101, 61, 153), QColor(153, 61, 75) });


/* ******************************************************************************************** *
 * Implementation of PlotWindow
 * ******************************************************************************************** */
PlotWindow::PlotWindow(QWidget *parent)
  : QMainWindow(parent), _plot(0)
{
  _plot = new QCustomPlot(this);
  _plot->setInteraction(QCP::iRangeZoom, true);
  _plot->setInteraction(QCP::iRangeDrag, true);
  setCentralWidget(_plot);

  QToolBar *toolbar = new QToolBar();
  toolbar->addAction(tr("Save"), this, SLOT(onSave()));
  this->addToolBar(toolbar);
}

void
PlotWindow::onSave() {
  QString filename = QFileDialog::getSaveFileName(this, tr("Save plot as ..."), "",
                                                  tr("Image Format (*.png *.pdf)"));
  if (filename.isEmpty())
    return;
  QFileInfo info(filename);
  if ("png" == info.suffix())
    _plot->savePng(filename);
  else if ("pdf" == info.suffix())
    _plot->savePdf(filename);
}


/* ******************************************************************************************** *
 * Implementation of MarginalPlotWindow
 * ******************************************************************************************** */
MarginalPlotWindow::MarginalPlotWindow(double tmin, double tmax, size_t nstep, const QVector<stochbb::Var> &vars, QWidget *parent)
  : PlotWindow(parent), _tmin(tmin), _tmax(tmax), _nstep(nstep), _vars(vars)
{
  double t = _tmin, dt = (_tmax-_tmin)/_nstep;
  Eigen::VectorXd F(_nstep), T(_nstep);
  for (size_t i=0; i<_nstep; i++, t+=dt)
    T(i) = t;
  double ymax = 0;
  for (size_t i=0; i<_vars.size(); i++) {
    _vars[i].density().eval(tmin, tmax, F);
    ymax = std::max(ymax, F.maxCoeff());
    QCPGraph *graph = _plot->addGraph();
    for (size_t j=0; j<nstep; j++)
      graph->addData(T[j], F[j]);
    // Set pen color;
    QPen pen = graph->pen();
    pen.setColor(colors[i % colors.size()]);
    pen.setWidth(2);
    graph->setPen(pen);
    if (_vars[i].name().size())
      graph->setName(QString::fromStdString(_vars[i].name()));
    graph->addToLegend();
  }
  _plot->legend->setVisible(true);

  _plot->xAxis->setRange(tmin, tmax);
  _plot->yAxis->setRange(0, ymax);
  _plot->replot();
}


/* ******************************************************************************************** *
 * Implementation of ScatterPlotWindow
 * ******************************************************************************************** */
ScatterPlotWindow::ScatterPlotWindow(size_t nsamples, const stochbb::Var &X, const stochbb::Var &Y, QWidget *parent)
  : PlotWindow(parent), _samples(nsamples,2)
{
  stochbb::ExactSampler sampler(X,Y);
  sampler.sample(_samples);

  if (X.name().size())
    _plot->xAxis->setLabel(QString::fromStdString(X.name()));
  if (Y.name().size())
    _plot->yAxis->setLabel(QString::fromStdString(Y.name()));

  QCPGraph *graph = _plot->addGraph();
  graph->setScatterStyle(QCPScatterStyle::ssPlus);
  graph->setLineStyle(QCPGraph::lsNone);

  double xmin = _samples.col(0).minCoeff(), xmax = _samples.col(0).maxCoeff();
  double ymin = _samples.col(1).minCoeff(), ymax = _samples.col(1).maxCoeff();

  for (int i=0; i<_samples.rows(); i++) {
    graph->addData(_samples(i, 0), _samples(i, 1));
  }

  _plot->xAxis->setRange(xmin, xmax);
  _plot->yAxis->setRange(ymin, ymax);
  _plot->replot();
}


/* ******************************************************************************************** *
 * Implementation of KDEPlotWindow
 * ******************************************************************************************** */
KDEPlotWindow::KDEPlotWindow(size_t nsamples, const QVector<stochbb::Var> &vars, QWidget *parent)
  : PlotWindow(parent), _nsamples(nsamples),  _vars(vars)
{
  Eigen::MatrixXd samples(_nsamples, _vars.size());
  stochbb::ExactSampler sampler(_vars.toStdVector());
  sampler.sample(samples);

  _densities.reserve(_vars.size());
  for (int i=0; i<samples.cols(); i++) {
    _densities.push_back(new KDE(samples.col(i)));
  }

  double min = _densities.first()->min();
  double max = _densities.first()->max();
  for (size_t i=1; i<_densities.size(); i++) {
    min = std::min(min, _densities[i]->min());
    max = std::max(max, _densities[i]->max());
  }

  size_t nstep = 200;
  double t = min, dt = (max-min)/nstep;
  Eigen::VectorXd F(nstep), T(nstep);
  for (size_t i=0; i<nstep; i++, t+=dt)
    T(i) = t;

  double ymax = 0;
  for (size_t i=0; i<_vars.size(); i++) {
    QCPGraph *graph = _plot->addGraph();
    for (size_t j=0; j<nstep; j++) {
      double y = _densities[i]->eval(T[j]);
      ymax = std::max(ymax, y);
      graph->addData(T[j], y);
    }
    // Set pen color;
    QPen pen = graph->pen();
    pen.setColor(colors[i % colors.size()]);
    pen.setWidth(2);
    graph->setPen(pen);
    if (_vars[i].name().size())
      graph->setName(QString::fromStdString(_vars[i].name()));
    graph->addToLegend();
  }
  _plot->legend->setVisible(true);
  _plot->xAxis->setRange(min, max);
  _plot->yAxis->setRange(0, ymax);
  _plot->replot();
}

KDEPlotWindow::~KDEPlotWindow() {
  for (size_t i=0; i<_densities.size(); i++) {
    delete _densities[i];
  }
  _densities.clear();
}


/* ******************************************************************************************** *
 * Implementation of KDE
 * ******************************************************************************************** */
KDE::KDE(const Eigen::Ref<Eigen::VectorXd> &samples)
  : _min(samples.minCoeff()), _max(samples.maxCoeff()), _bw(1), _samples(samples)
{
  double mu = _samples.mean();
  double sigma2 = (_samples.array()*_samples.array()).mean() - mu*mu;
  _bw = std::exp( std::log(1.06) + std::log(sigma2)/2 - std::log(_samples.size())/5 );
  _min -= 3*_bw;
  _max += 3*_bw;
}

double
KDE::eval(double x) const {
  double res = 0;
  for (int i=0; i<_samples.size(); i++)
    res += std::exp(-0.5*(x-_samples(i))*(x-_samples(i))/(2*_bw*_bw));
  return res / (std::sqrt(2*M_PI)*_samples.size()*_bw);
}

double
KDE::min() const {
  return _min;
}

double
KDE::max() const {
  return _max;
}
