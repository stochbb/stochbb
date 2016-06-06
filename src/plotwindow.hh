#ifndef PLOTWINDOW_HH
#define PLOTWINDOW_HH

#include <QMainWindow>
#include "lib/api.hh"
#include "qcustomplot.hh"

class PlotWindow: public QMainWindow
{
  Q_OBJECT

public:
  PlotWindow(QWidget *parent=0);

protected slots:
  void onSave();

protected:
  QCustomPlot *_plot;
};


class MarginalPlotWindow: public PlotWindow
{
  Q_OBJECT

public:
  MarginalPlotWindow(double tmin, double tmax, size_t nstep, const QVector<stochbb::Var> &vars,
                     QWidget *parent=0);

protected:
  double _tmin;
  double _tmax;
  size_t _nstep;
  QVector<stochbb::Var> _vars;
};


class ScatterPlotWindow: public PlotWindow
{
  Q_OBJECT

public:
  ScatterPlotWindow(size_t nsamples, const stochbb::Var &X, const stochbb::Var &Y, QWidget *parent=0);

protected:
  Eigen::MatrixXd _samples;
};


class KDE
{
public:
  explicit KDE(const Eigen::Ref<Eigen::VectorXd> &samples);

  double eval(double x) const;
  double min() const;
  double max() const;

protected:
  double _min, _max, _bw;
  Eigen::Ref<Eigen::VectorXd> _samples;
};


class KDEPlotWindow: public PlotWindow
{
  Q_OBJECT

public:
  KDEPlotWindow(size_t nsamples, const QVector<stochbb::Var> &vars,
                QWidget *parent=0);
  virtual ~KDEPlotWindow();

protected:
  size_t _nsamples;
  size_t _nbins;
  QVector<stochbb::Var> _vars;
  QVector<KDE *> _densities;
};

#endif // PLOTWINDOW_HH
