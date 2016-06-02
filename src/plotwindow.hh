#ifndef PLOTWINDOW_HH
#define PLOTWINDOW_HH

#include <QWidget>
#include "lib/api.hh"
#include "qcustomplot.hh"


class PlotWindow: public QWidget
{
  Q_OBJECT

public:
  PlotWindow(double tmin, double tmax, size_t nstep, const QVector<stochbb::Var> &vars,
             QWidget *parent=0);

protected:
  double _tmin;
  double _tmax;
  size_t _nstep;
  QVector<stochbb::Var> _vars;
};


#endif // PLOTWINDOW_HH
