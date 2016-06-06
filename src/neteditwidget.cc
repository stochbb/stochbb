#include "neteditwidget.hh"
#include "network.hh"
#include "nodes.hh"

#include <QVBoxLayout>
#include <QScrollArea>
#include <QFile>


NetEditWidget::NetEditWidget(QWidget *parent)
  : QWidget(parent), _netview(0)
{
  _netview = new Network();
  connect(_netview, SIGNAL(editNodeConfig(NodeBase*)), this, SLOT(onEditNodeConfig(NodeBase *)));

  QScrollArea *scroll = new QScrollArea();
  scroll->setWidget(_netview);
  scroll->setWidgetResizable(true);

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(scroll);
  layout->setSpacing(0);
  layout->setMargin(0);
  setLayout(layout);
}

Network *
NetEditWidget::network() {
  return _netview;
}

void
NetEditWidget::addDelay() {
  _netview->addNode(new DelayNode(_netview));
}

void
NetEditWidget::addRandomDelay() {
  _netview->addNode(new RandomDelayNode(_netview));
}

void
NetEditWidget::addGammaProc() {
  _netview->addNode(new GammaProcessNode(_netview));
}

void
NetEditWidget::addCompoundGammaProc() {
  _netview->addNode(new CompoundGammaProcessNode(_netview));
}

void
NetEditWidget::addInvGammaProc() {
  _netview->addNode(new InvGammaProcessNode(_netview));
}

void
NetEditWidget::addCompoundInvGammaProc() {
  _netview->addNode(new CompoundInvGammaProcessNode(_netview));
}

void
NetEditWidget::addWeibullProc() {
  _netview->addNode(new WeibullProcessNode(_netview));
}

void
NetEditWidget::addCompoundWeibullProc() {
  _netview->addNode(new CompoundWeibullProcessNode(_netview));
}

void
NetEditWidget::addMinimum() {
  _netview->addNode(new MinimumNode(_netview));
}

void
NetEditWidget::addMaximum() {
  _netview->addNode(new MaximumNode(_netview));
}

void
NetEditWidget::addInhibition() {
  InhibitionNode *inh = new InhibitionNode(_netview);
  JoinNode *join = new JoinNode(inh, _netview);
  _netview->addNode(inh);
  _netview->addNode(join);
}

void
NetEditWidget::addAffine() {
  _netview->addNode(new AffineNode(_netview));
}

void
NetEditWidget::addStimulus() {
  _netview->addNode(new TriggerNode(_netview));
}

void
NetEditWidget::addConstant() {
  _netview->addNode(new ConstantNode(_netview));
}

void
NetEditWidget::addUniformVar() {
  _netview->addNode(new UniformVarNode(_netview));
}

void
NetEditWidget::addNormalVar() {
  _netview->addNode(new NormalVarNode(_netview));
}

void
NetEditWidget::addCompoundNormalVar() {
  _netview->addNode(new CompoundNormalVarNode(_netview));
}

void
NetEditWidget::addGammaVar() {
  _netview->addNode(new GammaVarNode(_netview));
}

void
NetEditWidget::addCompoundGammaVar() {
  _netview->addNode(new CompoundGammaVarNode(_netview));
}

void
NetEditWidget::addInvGammaVar() {
  _netview->addNode(new InvGammaVarNode(_netview));
}

void
NetEditWidget::addCompoundInvGammaVar() {
  _netview->addNode(new CompoundInvGammaVarNode(_netview));
}

void
NetEditWidget::addWeibullVar() {
  _netview->addNode(new WeibullVarNode(_netview));
}

void
NetEditWidget::addCompoundWeibullVar() {
  _netview->addNode(new CompoundWeibullVarNode(_netview));
}

void
NetEditWidget::addMarginalPlot() {
  _netview->addNode(new MarginalPlotNode(_netview));
  _netview->update();
}

void
NetEditWidget::addScatterPlot() {
  _netview->addNode(new ScatterPlotNode(_netview));
}

void
NetEditWidget::addKDEPlot() {
  _netview->addNode(new KDEPlotNode(_netview));
}

void
NetEditWidget::removeSelected() {
  if (! _netview->selected())
    return;
  if (QNetNode *node = _netview->selectedNode()) {
    _netview->remNode(node);
  } else if (QNetEdge *edge = _netview->selectedEdge()) {
    _netview->remEdge(edge);
  }
}

void
NetEditWidget::zoomIn() {
  _netview->setScale(_netview->scale()+0.1);
}

void
NetEditWidget::zoomOut() {
  if (_netview->scale() > 0.1)
    _netview->setScale(_netview->scale()-0.1);
}

void
NetEditWidget::onEditNodeConfig(NodeBase *node) {
  if (QDialog::Accepted == NodeConfigDialog(node).exec())
    _netview->setModified(true);
}
