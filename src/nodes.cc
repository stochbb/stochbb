#include "nodes.hh"
#include "network.hh"
#include "assembler.hh"
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
  {"trigger",      (NodeBase *(*)(const QDomElement &node)) TriggerNode::fromXml},
  {"gammap",       (NodeBase *(*)(const QDomElement &node)) GammaProcessNode::fromXml},
  {"cgammap",      (NodeBase *(*)(const QDomElement &node)) CompoundGammaProcessNode::fromXml},
  {"invgammap",    (NodeBase *(*)(const QDomElement &node)) InvGammaProcessNode::fromXml},
  {"cinvgammap",   (NodeBase *(*)(const QDomElement &node)) CompoundInvGammaProcessNode::fromXml},
  {"weibullp",     (NodeBase *(*)(const QDomElement &node)) WeibullProcessNode::fromXml},
  {"cweibullp",    (NodeBase *(*)(const QDomElement &node)) CompoundWeibullProcessNode::fromXml},
  {"minimum",      (NodeBase *(*)(const QDomElement &node)) MinimumNode::fromXml},
  {"maximum",      (NodeBase *(*)(const QDomElement &node)) MaximumNode::fromXml},
  {"inhibition",   (NodeBase *(*)(const QDomElement &node)) InhibitionNode::fromXml},
  {"affine",       (NodeBase *(*)(const QDomElement &node)) AffineNode::fromXml},
  {"const",        (NodeBase *(*)(const QDomElement &node)) ConstantNode::fromXml},
  {"unifv",        (NodeBase *(*)(const QDomElement &node)) UniformVarNode::fromXml},
  {"normv",        (NodeBase *(*)(const QDomElement &node)) NormalVarNode::fromXml},
  {"cnormv",       (NodeBase *(*)(const QDomElement &node)) CompoundNormalVarNode::fromXml},
  {"gammav",       (NodeBase *(*)(const QDomElement &node)) GammaVarNode::fromXml},
  {"cgammav",      (NodeBase *(*)(const QDomElement &node)) CompoundGammaVarNode::fromXml},
  {"invgammav",    (NodeBase *(*)(const QDomElement &node)) InvGammaVarNode::fromXml},
  {"cinvgammav",   (NodeBase *(*)(const QDomElement &node)) CompoundInvGammaVarNode::fromXml},
  {"weibullv",     (NodeBase *(*)(const QDomElement &node)) WeibullVarNode::fromXml},
  {"cweibullv",    (NodeBase *(*)(const QDomElement &node)) CompoundWeibullVarNode::fromXml},
  {"marginalplot", (NodeBase *(*)(const QDomElement &node)) MarginalPlotNode::fromXml},
  {"scatterplot",  (NodeBase *(*)(const QDomElement &node)) ScatterPlotNode::fromXml},
  {"kdeplot",      (NodeBase *(*)(const QDomElement &node)) KDEPlotNode::fromXml}});


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

bool
NodeBase::needsPreprocessing(Assembler &assembler) const {
  return false;
}

bool
NodeBase::preprocess(Assembler &assembler) const {
  return false;
}

bool
NodeBase::processable(Assembler &assembler) const {
  // A node is processable if all its input sockets are connected to processed
  // output sockets (stored in _varTable).
  if (numSockets(QNetSocket::LEFT)) {
    // Iterate over all input sockets (on the left side)
    for (size_t i=0; i<numSockets(QNetSocket::LEFT); i++) {
      // get socket
      Socket *dest = dynamic_cast<Socket *>(socketAt(QNetSocket::LEFT, i));
      if ((! dest) || (! assembler.hasVariable(dest)))
        return false;
    }
  }
  return true;
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
 * Implementation of TriggerNode
 * ********************************************************************************************* */
TriggerNode::TriggerNode(Network *parent)
  : NodeBase("Trigger", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  _params.insert("time", 0);
}

QDomElement
TriggerNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","trigger");
  return node;
}

bool
TriggerNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  if (! out)
    return false;
  return assembler.addVariable(
        out, stochbb::delta(parameter("time"), label().toStdString()));
}

TriggerNode *
TriggerNode::fromXml(const QDomElement &node) {
  return  new TriggerNode();
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

bool
DelayNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  stochbb::Var in = assembler.sourceVar(this, "in");
  if (! out || in.isNull())
    return false;
  stochbb::Var res = stochbb::delta(parameter("delay"))+in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
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

bool
RandomDelayNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  stochbb::Var in = assembler.sourceVar(this, "in");
  stochbb::Var d  = assembler.sourceVar(this, "delay");
  if (! out || in.isNull() || d.isNull())
    return false;
  stochbb::Var res = d + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
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

bool
GammaProcessNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  stochbb::Var in = assembler.sourceVar(this, "in");
  if (!out || in.isNull())
    return false;
  stochbb::Var res = stochbb::gamma(parameter("k"), parameter("theta")) + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
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

bool
CompoundGammaProcessNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  stochbb::Var in = assembler.sourceVar(this, "in");
  stochbb::Var k = assembler.sourceVar(this, "k");
  stochbb::Var theta = assembler.sourceVar(this, "theta");
  if (!out || in.isNull() || k.isNull() || theta.isNull())
    return false;
  stochbb::Var res = stochbb::gamma(k, theta) + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
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

bool
InvGammaProcessNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  stochbb::Var in = assembler.sourceVar(this, "in");
  if (!out || in.isNull())
    return false;
  stochbb::Var res = stochbb::invgamma(parameter("alpha"), parameter("beta")) + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
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

bool
CompoundInvGammaProcessNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  stochbb::Var in = assembler.sourceVar(this, "in");
  stochbb::Var alpha = assembler.sourceVar(this, "alpha");
  stochbb::Var beta = assembler.sourceVar(this, "beta");
  if (!out || in.isNull() || alpha.isNull() || beta.isNull())
    return false;
  stochbb::Var res = stochbb::invgamma(alpha, beta) + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
}

CompoundInvGammaProcessNode *
CompoundInvGammaProcessNode::fromXml(const QDomElement &node) {
  return new CompoundInvGammaProcessNode();
}


/* ********************************************************************************************* *
 * Implementation of WeibullProcessNode
 * ********************************************************************************************* */
WeibullProcessNode::WeibullProcessNode(Network *parent)
  : NodeBase("Weibull", parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));

  this->_params.insert("k", 1);
  this->_params.insert("lambda", 1);
}

QDomElement
WeibullProcessNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","weibullp");
  return node;
}

bool
WeibullProcessNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  stochbb::Var in = assembler.sourceVar(this, "in");
  if (!out || in.isNull())
    return false;
  stochbb::Var res = stochbb::weibull(parameter("k"), parameter("lambda")) + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
}

WeibullProcessNode *
WeibullProcessNode::fromXml(const QDomElement &node) {
  return new WeibullProcessNode();
}


/* ********************************************************************************************* *
 * Implementation of CompoundWeibullProcessNode
 * ********************************************************************************************* */
CompoundWeibullProcessNode::CompoundWeibullProcessNode(Network *parent)
  : NodeBase("Weibull", parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "k", "k", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "lambda", QChar(0x03BB), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
}

QDomElement
CompoundWeibullProcessNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","cweibullp");
  return node;
}

bool
CompoundWeibullProcessNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  stochbb::Var in = assembler.sourceVar(this, "in");
  stochbb::Var k = assembler.sourceVar(this, "k");
  stochbb::Var lambda = assembler.sourceVar(this, "lambda");
  if (!out || in.isNull() || k.isNull() || lambda.isNull())
    return false;
  stochbb::Var res = stochbb::weibull(k, lambda) + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
}

CompoundWeibullProcessNode *
CompoundWeibullProcessNode::fromXml(const QDomElement &node) {
  return new CompoundWeibullProcessNode();
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

bool
MinimumNode::assemble(Assembler &assembler) const {
  stochbb::Var X = assembler.sourceVar(this, "X");
  stochbb::Var Y = assembler.sourceVar(this, "Y");
  Socket *out = assembler.socket(this, "out");
  if (!out || X.isNull() || Y.isNull())
    return false;
  stochbb::Var res = stochbb::minimum(X, Y);
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
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

bool
MaximumNode::assemble(Assembler &assembler) const {
  stochbb::Var X = assembler.sourceVar(this, "X");
  stochbb::Var Y = assembler.sourceVar(this, "Y");
  Socket *out = assembler.socket(this, "out");
  if (!out || X.isNull() || Y.isNull())
    return false;
  stochbb::Var res = stochbb::maximum(X, Y);
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
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

bool
InhibitionNode::needsPreprocessing(Assembler &assembler) const {
  return (! assembler.hasVariable(assembler.socket(this, "Aout"))) ||
      (! assembler.hasVariable(assembler.socket(this, "Bout")));
}

bool
InhibitionNode::preprocess(Assembler &assembler) const {
  if ((! assembler.socket(this, "Aout")) || (! assembler.socket(this, "Bout")))
    return false;
  if (! assembler.hasVariable(socket("Aout")))
    assembler.addVariable(socket("Aout"), stochbb::delta(0));
  if (! assembler.hasVariable(socket("Bout")))
    assembler.addVariable(socket("Bout"), stochbb::delta(0));
  return true;
}

bool
InhibitionNode::assemble(Assembler &assembler) const {
  stochbb::Var X = assembler.sourceVar(this, "inX");
  stochbb::Var Y = assembler.sourceVar(this, "inY");
  stochbb::Var A = assembler.sourceVar(this, "Ain");
  stochbb::Var B = assembler.sourceVar(this, "Bin");
  Socket *out = assembler.socket(this, "out");
  if ((!out) || X.isNull() || Y.isNull() || A.isNull() || B.isNull())
    return false;
  stochbb::Var res = stochbb::condchain(X, Y, A, B);
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
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

bool
AffineNode::assemble(Assembler &assembler) const {
  stochbb::Var in = assembler.sourceVar(this, "in");
  Socket *out = assembler.socket(this, "out");
  if (! out || in.isNull())
    return false;
  stochbb::Var res = parameter("scale") * in + parameter("shift");
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
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

bool
ConstantNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  if (! out)
    return false;
  return assembler.addVariable(out, stochbb::delta(parameter("value"), label().toStdString()));
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

bool
GammaVarNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  if (! out)
    return false;
  return assembler.addVariable(
        out, stochbb::gamma(parameter("k"), parameter("theta"), label().toStdString()));
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
  this->addSocket(new Socket(QNetSocket::LEFT, "k", tr("k"), this));
  this->addSocket(new Socket(QNetSocket::LEFT, "theta", QChar(0x0398), this));
}

QDomElement
CompoundGammaVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","cgammav");
  return node;
}

bool
CompoundGammaVarNode::assemble(Assembler &assembler) const {
  stochbb::Var k = assembler.sourceVar(this, "k");
  stochbb::Var theta = assembler.sourceVar(this, "theta");
  Socket *out = assembler.socket(this, "out");
  if (!out || k.isNull() || theta.isNull())
    return false;
  return assembler.addVariable(out, stochbb::gamma(k, theta, label().toStdString()));
}

CompoundGammaVarNode *
CompoundGammaVarNode::fromXml(const QDomElement &node) {
  return new CompoundGammaVarNode();
}


/* ********************************************************************************************* *
 * Implementation of InvGammaVarNode
 * ********************************************************************************************* */
InvGammaVarNode::InvGammaVarNode(Network *parent)
  : NodeBase(QString("1/")+QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("alpha", 1);
  this->_params.insert("beta", 1);
}

QDomElement
InvGammaVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","invgammav");
  return node;
}

bool
InvGammaVarNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  if (! out)
    return false;
  return assembler.addVariable(
        out, stochbb::invgamma(parameter("alpha"), parameter("beta"), label().toStdString()));
}

InvGammaVarNode *
InvGammaVarNode::fromXml(const QDomElement &node) {
  return new InvGammaVarNode();
}


/* ********************************************************************************************* *
 * Implementation of CompoundInvGammaVarNode
 * ********************************************************************************************* */
CompoundInvGammaVarNode::CompoundInvGammaVarNode(Network *parent)
  : NodeBase(QString("1/")+QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "alpha", QChar(0x03B1), this));
  this->addSocket(new Socket(QNetSocket::LEFT, "beta", QChar(0x03B2), this));
}

QDomElement
CompoundInvGammaVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","cinvgammav");
  return node;
}

bool
CompoundInvGammaVarNode::assemble(Assembler &assembler) const {
  stochbb::Var alpha = assembler.sourceVar(this, "alpha");
  stochbb::Var beta = assembler.sourceVar(this, "beta");
  Socket *out = assembler.socket(this, "out");
  if (!out || alpha.isNull() || beta.isNull())
    return false;
  return assembler.addVariable(out, stochbb::invgamma(alpha, beta, label().toStdString()));
}

CompoundInvGammaVarNode *
CompoundInvGammaVarNode::fromXml(const QDomElement &node) {
  return new CompoundInvGammaVarNode();
}


/* ********************************************************************************************* *
 * Implementation of WeibullVarNode
 * ********************************************************************************************* */
WeibullVarNode::WeibullVarNode(Network *parent)
  : NodeBase("Weibull", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("k", 1);
  this->_params.insert("lambda", 1);
}

QDomElement
WeibullVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","weibullv");
  return node;
}

bool
WeibullVarNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  if (! out)
    return false;
  return assembler.addVariable(
        out, stochbb::weibull(parameter("k"), parameter("lambda"), label().toStdString()));
}

WeibullVarNode *
WeibullVarNode::fromXml(const QDomElement &node) {
  return new WeibullVarNode();
}


/* ********************************************************************************************* *
 * Implementation of CompoundWeibullVarNode
 * ********************************************************************************************* */
CompoundWeibullVarNode::CompoundWeibullVarNode(Network *parent)
  : NodeBase("Weibull", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "k", "k", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "lambda", QChar(0x03BB), this));
}

QDomElement
CompoundWeibullVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","cweibullv");
  return node;
}

bool
CompoundWeibullVarNode::assemble(Assembler &assembler) const {
  stochbb::Var k = assembler.sourceVar(this, "k");
  stochbb::Var lambda = assembler.sourceVar(this, "lambda");
  Socket *out = assembler.socket(this, "out");
  if (!out || k.isNull() || lambda.isNull())
    return false;
  return assembler.addVariable(out, stochbb::weibull(k, lambda, label().toStdString()));
}

CompoundWeibullVarNode *
CompoundWeibullVarNode::fromXml(const QDomElement &node) {
  return new CompoundWeibullVarNode();
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

bool
UniformVarNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  if (! out)
    return false;
  return assembler.addVariable(
        out, stochbb::uniform(parameter("min"), parameter("max"), label().toStdString()));
}

UniformVarNode *
UniformVarNode::fromXml(const QDomElement &node) {
  return new UniformVarNode();
}


/* ********************************************************************************************* *
 * Implementation of NormalVarNode
 * ********************************************************************************************* */
NormalVarNode::NormalVarNode(Network *parent)
  : NodeBase("N", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("mu", 0);
  this->_params.insert("sigma", 1);
}

QDomElement
NormalVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","normv");
  return node;
}

bool
NormalVarNode::assemble(Assembler &assembler) const {
  Socket *out = assembler.socket(this, "out");
  if (! out)
    return false;
  return assembler.addVariable(
        out, stochbb::normal(parameter("mu"), parameter("sigma"), label().toStdString()));
}

NormalVarNode *
NormalVarNode::fromXml(const QDomElement &node) {
  return new NormalVarNode();
}


/* ********************************************************************************************* *
 * Implementation of CompoundNormalVarNode
 * ********************************************************************************************* */
CompoundNormalVarNode::CompoundNormalVarNode(Network *parent)
  : NodeBase("Weibull", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->addSocket(new Socket(QNetSocket::LEFT, "mu", QChar(0x03BC), this));
  this->addSocket(new Socket(QNetSocket::LEFT, "sigma", QChar(0x03C3), this));
}

QDomElement
CompoundNormalVarNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","cnormalv");
  return node;
}

bool
CompoundNormalVarNode::assemble(Assembler &assembler) const {
  stochbb::Var mu = assembler.sourceVar(this, "mu");
  stochbb::Var sigma = assembler.sourceVar(this, "sigma");
  Socket *out = assembler.socket(this, "out");
  if (!out || mu.isNull() || sigma.isNull())
    return false;
  return assembler.addVariable(out, stochbb::normal(mu, sigma, label().toStdString()));
}

CompoundNormalVarNode *
CompoundNormalVarNode::fromXml(const QDomElement &node) {
  return new CompoundNormalVarNode();
}


/* ********************************************************************************************* *
 * Implementation of OutputNode
 * ********************************************************************************************* */
OutputNode::OutputNode(const QString &label, QNetView *parent)
  : NodeBase(label, parent)
{
  // pass...
}

bool
OutputNode::assemble(Assembler &assembler) const {
  return true;
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

  MarginalPlotWindow *plot = new MarginalPlotWindow(tmin, tmax, nstep, vars);
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


/* ********************************************************************************************* *
 * Implementation of ScatterPlotNode
 * ********************************************************************************************* */
ScatterPlotNode::ScatterPlotNode(Network *parent)
  : OutputNode("Scatter Plot", parent)
{
  addSocket(new Socket(QNetSocket::LEFT, "X", "X", this));
  addSocket(new Socket(QNetSocket::LEFT, "Y", "Y", this));
  _params.insert("samples", 1000);
}

void
ScatterPlotNode::execute(const QHash<Socket *, stochbb::Var> &vartable) {
  stochbb::Var X = vartable[socket("X")];
  stochbb::Var Y = vartable[socket("Y")];

  size_t samples = parameter("samples") > 0 ? parameter("samples") : 1000;

  ScatterPlotWindow *plot = new ScatterPlotWindow(samples, X, Y);
  plot->resize(480, 320);
  plot->setWindowTitle(this->label());
  plot->show();
}

QDomElement
ScatterPlotNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type", "scatterplot");
  return node;
}

ScatterPlotNode *
ScatterPlotNode::fromXml(const QDomElement &node) {
  return new ScatterPlotNode();
}


/* ********************************************************************************************* *
 * Implementation of KDEPlotNode
 * ********************************************************************************************* */
KDEPlotNode::KDEPlotNode(Network *parent)
  : OutputNode("KDE Plot", parent)
{
  _params.insert("graphs", 0);
  _params.insert("samples", 100);
}

bool
KDEPlotNode::setParameter(const QString &name, double value) {
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
KDEPlotNode::execute(const QHash<Socket *, stochbb::Var> &vartable) {
  if (0 == numSockets(QNetSocket::LEFT))
    return;
  QVector<stochbb::Var> vars;
  for (size_t i=0; i<numSockets(QNetSocket::LEFT); i++) {
    stochbb::Var X = vartable[socket(QString::number(i+1))];
    if (X.isNull())
      continue;
    vars.push_back(X);
  }

  size_t nsample = parameter("samples") > 0 ? parameter("samples") : 100;

  KDEPlotWindow *plot = new KDEPlotWindow(nsample, vars);
  plot->resize(480, 320);
  plot->setWindowTitle(this->label());
  plot->show();
}

QDomElement
KDEPlotNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type", "kdeplot");
  return node;
}

KDEPlotNode *
KDEPlotNode::fromXml(const QDomElement &node) {
  return new KDEPlotNode();
}


