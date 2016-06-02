#include "nodes.hh"
#include "network.hh"
#include "plotwindow.hh"
#include <sstream>
#include <QFormLayout>
#include <QLineEdit>
#include <QDoubleValidator>
#include <QLabel>
#include <QDialogButtonBox>


/* ********************************************************************************************* *
 * Implementation of NodeConfigDialog
 * ********************************************************************************************* */
NodeConfigDialog::NodeConfigDialog(NodeBase *node)
  : QDialog(), _node(node)
{
  QFormLayout *nodeprops = new QFormLayout();
  _label = new QLineEdit(_node->label());
  nodeprops->addRow(tr("Label"), _label);

  QFormLayout *param_layout = 0;
  if (node->hasParameters())
    param_layout = new QFormLayout();

  const QHash<QString, double> &params = node->parameters();
  QHash<QString, double>::const_iterator param = params.begin();
  for (; param != params.end(); param++) {
    QLineEdit *p = new QLineEdit(QString::number(param.value()));
    p->setValidator(new QDoubleValidator());
    _params.insert(param.key(), p);
    param_layout->addRow(param.key(), p);
  }

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addLayout(nodeprops);
  layout->addWidget(new QLabel(tr("Parameters:")));
  if (param_layout)
    layout->addLayout(param_layout);
  else
    layout->addWidget(new QLabel(tr("This node has no parameters.")));

  QDialogButtonBox *bb = new QDialogButtonBox(QDialogButtonBox::Cancel | QDialogButtonBox::Ok);
  connect(bb, SIGNAL(accepted()), this, SLOT(apply()));
  connect(bb, SIGNAL(rejected()), this, SLOT(reject()));
  layout->addWidget(bb);

  setLayout(layout);
}

void
NodeConfigDialog::apply() {
  _node->setLabel(_label->text());
  QHash<QString, QLineEdit *>::iterator item = _params.begin();
  for (; item != _params.end(); item++) {
    _node->setParameter(item.key(), item.value()->text().toDouble());
  }
  accept();
}


/* ********************************************************************************************* *
 * Implementation of Socket
 * ********************************************************************************************* */
Socket::Socket(Side side, const QString &name, const QString &label, QNetNode *node)
  : QNetSocket(side, label, node), _name(name)
{
  // pass...
}

const QString &
Socket::name() const {
  return _name;
}


/* ********************************************************************************************* *
 * Implementation of NodeBase
 * ********************************************************************************************* */
QHash<QString, NodeBase *(*)(const QDomElement &node)>
NodeBase::_factoryFunctions({
  {"delay",        (NodeBase *(*)(const QDomElement &node)) DelayNode::fromXml},
  {"rdelay",       (NodeBase *(*)(const QDomElement &node)) RandomDelayNode::fromXml},
  {"stimulus",     (NodeBase *(*)(const QDomElement &node)) StimulusNode::fromXml},
  {"gammap",       (NodeBase *(*)(const QDomElement &node)) GammaProcessNode::fromXml},
  {"cgammap",      (NodeBase *(*)(const QDomElement &node)) CompoundGammaProcessNode::fromXml},
  {"invgammap",    (NodeBase *(*)(const QDomElement &node)) InvGammaProcessNode::fromXml},
  {"cinvgammap",   (NodeBase *(*)(const QDomElement &node)) CompoundInvGammaProcessNode::fromXml},
  {"minimum",      (NodeBase *(*)(const QDomElement &node)) MinimumNode::fromXml},
  {"maximum",      (NodeBase *(*)(const QDomElement &node)) MaximumNode::fromXml},
  {"inhibition",   (NodeBase *(*)(const QDomElement &node)) InhibitionNode::fromXml},
  {"affine",       (NodeBase *(*)(const QDomElement &node)) AffineNode::fromXml},
  {"const",        (NodeBase *(*)(const QDomElement &node)) ConstantNode::fromXml},
  {"gammav",       (NodeBase *(*)(const QDomElement &node)) GammaVarNode::fromXml},
  {"cgammav",      (NodeBase *(*)(const QDomElement &node)) CompoundGammaVarNode::fromXml},
  {"unifv",        (NodeBase *(*)(const QDomElement &node)) UniformVarNode::fromXml},
  {"marginalplot", (NodeBase *(*)(const QDomElement &node)) MarginalPlotNode::fromXml} });


NodeBase::NodeBase(const QString &label, QNetView *parent)
  : QNetNode(label, parent)
{
  std::stringstream buffer; buffer << this;
  _id = QString(buffer.str().c_str());
}

const QString &
NodeBase::id() const {
  return _id;
}

bool
NodeBase::hasParameters() const {
  return 0 != _params.size();
}

bool
NodeBase::hasParameter(const QString &name) const {
  return _params.contains(name);
}

double
NodeBase::parameter(const QString &name) const {
  return _params[name];
}

const QHash<QString, double> &
NodeBase::parameters() const {
  return _params;
}

bool
NodeBase::setParameter(const QString &name, double value) {
  if (! _params.contains(name))
    return false;
  _params[name] = value;
  return true;
}

bool
NodeBase::hasSocket(const QString &name) const {
  return _sockets.contains(name);
}

Socket *
NodeBase::socket(const QString &name) const {
  return _sockets[name];
}

void
NodeBase::addSocket(QNetSocket *obj) {
  QNetNode::addSocket(obj);
  if (Socket *socket = dynamic_cast<Socket *>(obj)) {
    _sockets.insert(socket->name(), socket);
  }
}

QDomElement
NodeBase::serialize(QDomDocument &doc) const {
  QDomElement node = doc.createElement("node");
  node.setAttribute("id", _id);
  node.setAttribute("x", position().x());
  node.setAttribute("y", position().y());
  node.setAttribute("label", label());

  QHash<QString, double>::const_iterator param = _params.begin();
  for (; param != _params.end(); param++) {
    QDomElement pnode = doc.createElement("parameter");
    pnode.setAttribute("name", param.key());
    pnode.appendChild(doc.createTextNode(QString::number(param.value())));
    node.appendChild(pnode);
  }

  return node;
}

NodeBase *
NodeBase::fromXml(const QDomElement &node) {
  // Check if node has 'type' attribute
  if (! node.hasAttribute("type"))
    return 0;

  // Dispatch by type
  NodeBase *obj = 0;
  if (_factoryFunctions.contains(node.attribute("type"))) {
    obj = _factoryFunctions[node.attribute("type")](node);
  }

  if (obj) {
    QPoint pos(obj->position());
    if (node.hasAttribute("x"))
      pos.setX(node.attribute("x").toInt());
    if (node.hasAttribute("y"))
      pos.setY(node.attribute("y").toInt());
    obj->setPosition(pos);
    if (node.hasAttribute("label"))
      obj->setLabel(node.attribute("label"));

    QDomElement param = node.firstChildElement("parameter");
    for (; ! param.isNull(); param = param.nextSiblingElement("parameter")) {
      if (param.hasAttribute("name"))
        obj->setParameter(param.attribute("name"), param.text().toDouble());
    }
  }

  return obj;
}


/* ********************************************************************************************* *
 * Implementation of StimulusNode
 * ********************************************************************************************* */
StimulusNode::StimulusNode(Network *parent)
  : NodeBase("Stimulus", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  _params.insert("time", 0);
}

QDomElement
StimulusNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","stimulus");
  return node;
}

StimulusNode *
StimulusNode::fromXml(const QDomElement &node) {
  return  new StimulusNode();
}



/* ********************************************************************************************* *
 * Implementation of DelayNode
 * ********************************************************************************************* */
DelayNode::DelayNode(Network *parent)
  : NodeBase(QChar(0x03B4), parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));

  this->_params.insert("delay", 0);
}

QDomElement
DelayNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","delay");
  return node;
}

DelayNode *
DelayNode::fromXml(const QDomElement &node) {
  return new DelayNode();
}


/* ********************************************************************************************* *
 * Implementation of RandomDelayNode
 * ********************************************************************************************* */
RandomDelayNode::RandomDelayNode(Network *parent)
  : NodeBase(QChar(0x03B4), parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "delay", QChar(0x0394), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
}

QDomElement
RandomDelayNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","rdelay");
  return node;
}

RandomDelayNode *
RandomDelayNode::fromXml(const QDomElement &node) {
  return new RandomDelayNode();
}


/* ********************************************************************************************* *
 * Implementation of GammaProcessNode
 * ********************************************************************************************* */
GammaProcessNode::GammaProcessNode(Network *parent)
  : NodeBase(QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));

  this->_params.insert("k", 1);
  this->_params.insert("theta", 1);
}

QDomElement
GammaProcessNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","gammap");
  return node;
}

GammaProcessNode *
GammaProcessNode::fromXml(const QDomElement &node) {
  return new GammaProcessNode();
}


/* ********************************************************************************************* *
 * Implementation of GammaProcessNode
 * ********************************************************************************************* */
CompoundGammaProcessNode::CompoundGammaProcessNode(Network *parent)
  : NodeBase(QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "k", "k", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "theta", QChar(0x0398), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
}

QDomElement
CompoundGammaProcessNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","cgammap");
  return node;
}

CompoundGammaProcessNode *
CompoundGammaProcessNode::fromXml(const QDomElement &node) {
  return new CompoundGammaProcessNode();
}


/* ********************************************************************************************* *
 * Implementation of InvGammaProcessNode
 * ********************************************************************************************* */
InvGammaProcessNode::InvGammaProcessNode(Network *parent)
  : NodeBase(QString("1/")+QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));

  this->_params.insert("alpha", 1);
  this->_params.insert("beta", 1);
}

QDomElement
InvGammaProcessNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","invgammap");
  return node;
}

InvGammaProcessNode *
InvGammaProcessNode::fromXml(const QDomElement &node) {
  return new InvGammaProcessNode();
}


/* ********************************************************************************************* *
 * Implementation of CompoundInvGammaProcessNode
 * ********************************************************************************************* */
CompoundInvGammaProcessNode::CompoundInvGammaProcessNode(Network *parent)
  : NodeBase(QString("1/")+QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "alpha", QChar(0x03B1), this));
  this->addSocket(new Socket(QNetSocket::LEFT, "beta", QChar(0x03B2), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
}

QDomElement
CompoundInvGammaProcessNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","cinvgammap");
  return node;
}

CompoundInvGammaProcessNode *
CompoundInvGammaProcessNode::fromXml(const QDomElement &node) {
  return new CompoundInvGammaProcessNode();
}


/* ********************************************************************************************* *
 * Implementation of MinimumNode
 * ********************************************************************************************* */
MinimumNode::MinimumNode(Network *parent)
  : NodeBase("Min", parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "X", tr("X"), this));
  this->addSocket(new Socket(QNetSocket::LEFT, "Y", tr("Y"), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
}

QDomElement
MinimumNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","minimum");
  return node;
}

MinimumNode *
MinimumNode::fromXml(const QDomElement &node) {
  return new MinimumNode();
}


/* ********************************************************************************************* *
 * Implementation of MaximumNode
 * ********************************************************************************************* */
MaximumNode::MaximumNode(Network *parent)
  : NodeBase("Max", parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "X", tr("X"), this));
  this->addSocket(new Socket(QNetSocket::LEFT, "Y", tr("Y"), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
}

QDomElement
MaximumNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","maximum");
  return node;
}

MaximumNode *
MaximumNode::fromXml(const QDomElement &node) {
  return new MaximumNode();
}


/* ********************************************************************************************* *
 * Implementation of InhibitionNode
 * ********************************************************************************************* */
InhibitionNode::InhibitionNode(Network *parent)
  : NodeBase("Inh", parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT,  "Ain", tr(">A"), this));
  this->addSocket(new Socket(QNetSocket::LEFT,  "inX", tr("X"), this));
  this->addSocket(new Socket(QNetSocket::LEFT,  "inY", tr("Y"), this));
  this->addSocket(new Socket(QNetSocket::LEFT,  "Bin", tr(">B"), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "Aout", tr("A>"), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", tr("out"), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "Bout", tr("B>"), this));
}

QDomElement
InhibitionNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","inhibition");
  return node;
}

InhibitionNode *
InhibitionNode::fromXml(const QDomElement &node) {
  return new InhibitionNode();
}


/* ********************************************************************************************* *
 * Implementation of AffineNode
 * ********************************************************************************************* */
AffineNode::AffineNode(Network *parent)
  : NodeBase("a X + b", parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("scale", 1);
  this->_params.insert("shift", 0);
}

QDomElement
AffineNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","affine");
  return node;
}

AffineNode *
AffineNode::fromXml(const QDomElement &node) {
  return new AffineNode();
}


/* ********************************************************************************************* *
 * Implementation of ConstantNode
 * ********************************************************************************************* */
ConstantNode::ConstantNode(Network *parent)
  : NodeBase("C", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("value", 0);
}

QDomElement
ConstantNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","const");
  return node;
}

ConstantNode *
ConstantNode::fromXml(const QDomElement &node) {
  return new ConstantNode();
}


/* ********************************************************************************************* *
 * Implementation of GammaVarNode
 * ********************************************************************************************* */
GammaVarNode::GammaVarNode(Network *parent)
  : NodeBase(QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("k", 1);
  this->_params.insert("theta", 1);
}

QDomElement
GammaVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","gammav");
  return node;
}

GammaVarNode *
GammaVarNode::fromXml(const QDomElement &node) {
  return new GammaVarNode();
}


/* ********************************************************************************************* *
 * Implementation of CompoundGammaVarNode
 * ********************************************************************************************* */
CompoundGammaVarNode::CompoundGammaVarNode(Network *parent)
  : NodeBase(QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "theta", tr("k"), this));
  this->addSocket(new Socket(QNetSocket::LEFT, "theta", QChar(0x0398), this));
}

QDomElement
CompoundGammaVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","cgammav");
  return node;
}

CompoundGammaVarNode *
CompoundGammaVarNode::fromXml(const QDomElement &node) {
  return new CompoundGammaVarNode();
}


/* ********************************************************************************************* *
 * Implementation of UniformVarNode
 * ********************************************************************************************* */
UniformVarNode::UniformVarNode(Network *parent)
  : NodeBase("U", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("min", 0);
  this->_params.insert("max", 1);
}

QDomElement
UniformVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","unifv");
  return node;
}

UniformVarNode *
UniformVarNode::fromXml(const QDomElement &node) {
  return new UniformVarNode();
}


/* ********************************************************************************************* *
 * Implementation of OutputNode
 * ********************************************************************************************* */
OutputNode::OutputNode(const QString &label, QNetView *parent)
  : NodeBase(label, parent)
{
  // pass...
}


/* ********************************************************************************************* *
 * Implementation of MarginalPlotNode
 * ********************************************************************************************* */
MarginalPlotNode::MarginalPlotNode(Network *parent)
  : OutputNode("Marginal Plot", parent)
{
  _params.insert("graphs", 0);
  _params.insert("min", 0);
  _params.insert("max", 1);
  _params.insert("steps", 100);
}

bool
MarginalPlotNode::setParameter(const QString &name, double value) {
  if ("graphs" == name) {
    if (value < numSockets(QNetSocket::LEFT))
      return false;
    size_t n = numSockets(QNetSocket::LEFT)+1;
    while (value > numSockets(QNetSocket::LEFT)) {
      addSocket(new Socket(QNetSocket::LEFT, QString::number(n), QString::number(n), this));
      n++;
    }
  }
  return NodeBase::setParameter(name, value);
}

void
MarginalPlotNode::execute(const QHash<Socket *, stochbb::Var> &vartable) {
  if (0 == numSockets(QNetSocket::LEFT))
    return;
  QVector<stochbb::Var> vars;
  for (size_t i=0; i<numSockets(QNetSocket::LEFT); i++) {
    stochbb::Var X = vartable[socket(QString::number(i+1))];
    if (X.isNull())
      continue;
    vars.push_back(X);
  }

  size_t nstep = parameter("steps") > 0 ? parameter("steps") : 100;
  double tmin = parameter("min"), tmax = parameter("max");

  PlotWindow *plot = new PlotWindow(tmin, tmax, nstep, vars);
  plot->resize(480, 320);
  plot->setWindowTitle(this->label());
  plot->show();
}

QDomElement
MarginalPlotNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type", "marginalplot");
  return node;
}

MarginalPlotNode *
MarginalPlotNode::fromXml(const QDomElement &node) {
  return new MarginalPlotNode();
}
