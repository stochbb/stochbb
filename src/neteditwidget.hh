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
  void addMinimum();
  void addMaximum();
  void addInhibition();
  void addAffine();
  void addStimulus();
  void addConstant();
  void addGammaVar();
  void addGammaCompoundVar();
  void addMarginalPlot();
  void removeSelected();
  void zoomIn();
  void zoomOut();

protected slots:
  void onEditNodeConfig(NodeBase *node);

protected:
  Network *_netview;
};

#endif // NETEDITWIDGET_HH
