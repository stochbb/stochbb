#ifndef NETEDITWIDGET_HH
#define NETEDITWIDGET_HH

#include <QWidget>

class Network;
class NodeBase;


class NetEditWidget : public QWidget
{
  Q_OBJECT

public:
  explicit NetEditWidget(QWidget *parent = 0);

  Network *network();

public slots:
  void addDelay();
  void addRandomDelay();
  void addGammaProc();
  void addCompoundGammaProc();
  void addInvGammaProc();
  void addCompoundInvGammaProc();
  void addWeibullProc();
  void addCompoundWeibullProc();
  void addMinimum();
  void addMaximum();
  void addInhibition();
  void addAffine();
  void addStimulus();
  void addConstant();
  void addUniformVar();
  void addNormalVar();
  void addCompoundNormalVar();
  void addGammaVar();
  void addCompoundGammaVar();
  void addInvGammaVar();
  void addCompoundInvGammaVar();
  void addWeibullVar();
  void addCompoundWeibullVar();
  void addMarginalPlot();
  void addScatterPlot();
  void addKDEPlot();
  void removeSelected();
  void zoomIn();
  void zoomOut();

protected slots:
  void onEditNodeConfig(NodeBase *node);

protected:
  Network *_netview;
};

#endif // NETEDITWIDGET_HH
