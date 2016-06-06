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
QHash<QString, NodeBase::nodeFactoryFunction>
NodeBase::_factoryFunctions({
  {"delay",        (NodeBase::nodeFactoryFunction) DelayNode::fromXml},
  {"rdelay",       (NodeBase::nodeFactoryFunction) RandomDelayNode::fromXml},
  {"trigger",      (NodeBase::nodeFactoryFunction) TriggerNode::fromXml},
  {"gammap",       (NodeBase::nodeFactoryFunction) GammaProcessNode::fromXml},
  {"cgammap",      (NodeBase::nodeFactoryFunction) CompoundGammaProcessNode::fromXml},
  {"invgammap",    (NodeBase::nodeFactoryFunction) InvGammaProcessNode::fromXml},
  {"cinvgammap",   (NodeBase::nodeFactoryFunction) CompoundInvGammaProcessNode::fromXml},
  {"weibullp",     (NodeBase::nodeFactoryFunction) WeibullProcessNode::fromXml},
  {"cweibullp",    (NodeBase::nodeFactoryFunction) CompoundWeibullProcessNode::fromXml},
  {"minimum",      (NodeBase::nodeFactoryFunction) MinimumNode::fromXml},
  {"maximum",      (NodeBase::nodeFactoryFunction) MaximumNode::fromXml},
  {"inhibition",   (NodeBase::nodeFactoryFunction) InhibitionNode::fromXml},
  {"join",         (NodeBase::nodeFactoryFunction) JoinNode::fromXml},
  {"affine",       (NodeBase::nodeFactoryFunction) AffineNode::fromXml},
  {"const",        (NodeBase::nodeFactoryFunction) ConstantNode::fromXml},
  {"unifv",        (NodeBase::nodeFactoryFunction) UniformVarNode::fromXml},
  {"normv",        (NodeBase::nodeFactoryFunction) NormalVarNode::fromXml},
  {"cnormv",       (NodeBase::nodeFactoryFunction) CompoundNormalVarNode::fromXml},
  {"gammav",       (NodeBase::nodeFactoryFunction) GammaVarNode::fromXml},
  {"cgammav",      (NodeBase::nodeFactoryFunction) CompoundGammaVarNode::fromXml},
  {"invgammav",    (NodeBase::nodeFactoryFunction) InvGammaVarNode::fromXml},
  {"cinvgammav",   (NodeBase::nodeFactoryFunction) CompoundInvGammaVarNode::fromXml},
  {"weibullv",     (NodeBase::nodeFactoryFunction) WeibullVarNode::fromXml},
  {"cweibullv",    (NodeBase::nodeFactoryFunction) CompoundWeibullVarNode::fromXml},
  {"marginalplot", (NodeBase::nodeFactoryFunction) MarginalPlotNode::fromXml},
  {"scatterplot",  (NodeBase::nodeFactoryFunction) ScatterPlotNode::fromXml},
  {"kdeplot",      (NodeBase::nodeFactoryFunction) KDEPlotNode::fromXml}});


NodeBase::NodeBase(const QString &label, QNetView *parent)
  : QNetNode(label, parent)
{
  std::stringstream buffer; buffer << this;
  _id = QString(buffer.str().c_str());
  _type = "<NodeBase::_type attribute unset>";
}

const QString &
NodeBase::id() const {
  return _id;
}

const QString &
NodeBase::type() const {
  return _type;
}

bool
NodeBase::hasParameters() const {
  return 0 != _params.size();
}

bool
NodeBase::hasParameter(const QString &name) const {
  return _params.contains(name);
}

Parameter
NodeBase::parameter(const QString &name) const {
  return _params[name];
}

const QHash<QString, Parameter> &
NodeBase::parameters() const {
  return _params;
}

void
NodeBase::addParameter(const QString &name, const Parameter &param) {
  _params[name] = param;
}

bool
NodeBase::setParameter(const QString &name, const Parameter &param) {
  if (! _params.contains(name))
    return false;
  if (_params[name].type() != param.type())
    return false;
  _params[name] = param;
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

  QHash<QString, Parameter>::const_iterator param = _params.begin();
  for (; param != _params.end(); param++) {
    QDomElement pnode = param.value().serialize(doc);
    pnode.setAttribute("name", param.key());
    node.appendChild(pnode);
  }

  if (hasDescription()) {
    QDomElement descr = doc.createElement("description");
    descr.appendChild(doc.createTextNode(description()));
    node.appendChild(descr);
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
NodeBase::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  // Check if node has 'type' attribute
  if (! node.hasAttribute("type")) {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Node " << node.tagName()
        << " has no 'type' attribute.";
    info.addError(text);
    return 0;
  }


  // Dispatch by type
  NodeBase *obj = 0;
  if (_factoryFunctions.contains(node.attribute("type"))) {
    obj = _factoryFunctions[node.attribute("type")](node, info, nodeTable);
  } else {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Unknown node type '"
        << node.attribute("type") << "'.";
    info.addError(text);
    return 0;
  }

  QPoint pos(obj->position());
  if (node.hasAttribute("x"))
    pos.setX(node.attribute("x").toInt());
  else {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Node has no 'x' attribute.";
    info.addWarning(text);
  }
  if (node.hasAttribute("y"))
    pos.setY(node.attribute("y").toInt());
  else {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Node has no 'y' attribute.";
    info.addWarning(text);
  }
  obj->setPosition(pos);
  if (node.hasAttribute("label"))
    obj->setLabel(node.attribute("label"));
  QDomElement descr = node.firstChildElement("description");
  if (descr.isElement())
    obj->setDescription(descr.text());

  QDomElement param = node.firstChildElement("parameter");
  for (; ! param.isNull(); param = param.nextSiblingElement("parameter")) {
    if (param.hasAttribute("name")) {
       if(! obj->setParameter(param.attribute("name"), Parameter::fromXml(param, info))) {
         QString text; QTextStream msg(&text);
         msg << "@line " << node.lineNumber() << ": Cannot set parameter '"
             << param.attribute("name") << "' to '"
             << param.text() << "'.";
         info.addWarning(text);
       }
    } else {
      QString text; QTextStream msg(&text);
      msg << "@line " << node.lineNumber() << ": 'parameter' element has no name.";
      info.addWarning(text);
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
  _params.insert("time", Parameter(0.0));
  _type = "trigger";
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
        out, stochbb::delta(parameter("time").asFloat(), label().toStdString()));
}

TriggerNode *
TriggerNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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

  this->_params.insert("delay", Parameter(0.0));

  _type = "fixed delay";
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
  stochbb::Var res = stochbb::delta(parameter("delay").asFloat())+in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
}

DelayNode *
DelayNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _type = "random delay";
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
RandomDelayNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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

  this->_params.insert("k", Parameter(1.0));
  this->_params.insert("theta", Parameter(1.0));

  _type = "gamma process";
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
  stochbb::Var res = stochbb::gamma(parameter("k").asFloat(), parameter("theta").asFloat()) + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
}

GammaProcessNode *
GammaProcessNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _type = "compound gamma process";
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
CompoundGammaProcessNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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

  this->_params.insert("alpha", Parameter(1.0));
  this->_params.insert("beta", Parameter(1.0));

  _type = "inverse gamma process";
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
  stochbb::Var res = stochbb::invgamma(parameter("alpha").asFloat(), parameter("beta").asFloat()) + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
}

InvGammaProcessNode *
InvGammaProcessNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _type = "compound inverse gamma process";
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
CompoundInvGammaProcessNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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

  this->_params.insert("k", Parameter(1.0));
  this->_params.insert("lambda", Parameter(1.0));

  _type = "weibull process";
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
  stochbb::Var res = stochbb::weibull(parameter("k").asFloat(), parameter("lambda").asFloat()) + in;
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
}

WeibullProcessNode *
WeibullProcessNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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

  _type = "compound weibull process";
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
CompoundWeibullProcessNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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

  _type = "minimum";
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
MinimumNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _type = "maximum";
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
MaximumNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new MaximumNode();
}


/* ********************************************************************************************* *
 * Implementation of InhibitionNode
 * ********************************************************************************************* */
InhibitionNode::InhibitionNode(Network *parent)
  : NodeBase("Inh", parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT,  "X", tr("X"), this));
  this->addSocket(new Socket(QNetSocket::LEFT,  "Y", tr("Y"), this));
  this->addSocket(new Socket(QNetSocket::RIGHT,  "Xout", tr("X"), this));
  this->addSocket(new Socket(QNetSocket::RIGHT,  "Yout", tr("Y"), this));
  _type = "inhibition";
}

QDomElement
InhibitionNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","inhibition");
  return node;
}

bool
InhibitionNode::assemble(Assembler &assembler) const {
  stochbb::Var X = assembler.sourceVar(this, "X");
  stochbb::Var Y = assembler.sourceVar(this, "Y");
  Socket *Xout = assembler.socket(this, "Xout");
  Socket *Yout = assembler.socket(this, "Yout");
  if ((!Xout) || (!Yout) || X.isNull() || Y.isNull())
    return false;
  assembler.addVariable(Xout, stochbb::delta(0));
  assembler.addVariable(Yout, stochbb::delta(0));
  return true;
}

InhibitionNode *
InhibitionNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new InhibitionNode();
}


/* ********************************************************************************************* *
 * Implementation of JoinNode
 * ********************************************************************************************* */
JoinNode::JoinNode(InhibitionNode *sibling, Network *parent)
  : NodeBase("Join", parent), _sibling(sibling)
{
  this->addSocket(new Socket(QNetSocket::LEFT,  "X", tr("X"), this));
  this->addSocket(new Socket(QNetSocket::LEFT,  "Y", tr("Y"), this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", tr("out"), this));
  _type = "join";

  connect(_sibling, SIGNAL(destroyed(QObject*)), this, SLOT(deleteLater()));
  connect(this, SIGNAL(destroyed(QObject*)), _sibling, SLOT(deleteLater()));
}

QDomElement
JoinNode::serialize(QDomDocument &doc) const {
  QDomElement node = NodeBase::serialize(doc);
  node.setAttribute("type","join");
  node.setAttribute("sibling", _sibling->id());
  return node;
}

bool
JoinNode::assemble(Assembler &assembler) const {
  stochbb::Var X = assembler.sourceVar(_sibling, "X");
  stochbb::Var Y = assembler.sourceVar(_sibling, "Y");
  stochbb::Var A = assembler.sourceVar(this, "X");
  stochbb::Var B = assembler.sourceVar(this, "Y");
  Socket *out = assembler.socket(this, "out");
  if ((!out) || X.isNull() || Y.isNull() || A.isNull() || B.isNull())
    return false;
  stochbb::Var res = stochbb::condchain(X, Y, A, B);
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
}

JoinNode *
JoinNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  if (! node.hasAttribute("sibling"))
    return 0;
  if (! nodeTable.contains(node.attribute("sibling")))
    return 0;
  InhibitionNode *inh = dynamic_cast<InhibitionNode *>(nodeTable[node.attribute("sibling")]);
  if (! inh)
    return 0;
  return new JoinNode(inh);
}


/* ********************************************************************************************* *
 * Implementation of AffineNode
 * ********************************************************************************************* */
AffineNode::AffineNode(Network *parent)
  : NodeBase("a X + b", parent)
{
  this->addSocket(new Socket(QNetSocket::LEFT, "in", "", this));
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("scale", Parameter(1.0));
  this->_params.insert("shift", Parameter(1.0));
  _type = "affine transform";
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
  stochbb::Var res = parameter("scale").asFloat() * in + parameter("shift").asFloat();
  res.setName(label().toStdString());
  return assembler.addVariable(out, res);
}

AffineNode *
AffineNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new AffineNode();
}


/* ********************************************************************************************* *
 * Implementation of ConstantNode
 * ********************************************************************************************* */
ConstantNode::ConstantNode(Network *parent)
  : NodeBase("C", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("value", Parameter(1.0));
  _type = "constant";
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
  return assembler.addVariable(out, stochbb::delta(parameter("value").asFloat(), label().toStdString()));
}

ConstantNode *
ConstantNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new ConstantNode();
}


/* ********************************************************************************************* *
 * Implementation of GammaVarNode
 * ********************************************************************************************* */
GammaVarNode::GammaVarNode(Network *parent)
  : NodeBase(QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("k", Parameter(1.0));
  this->_params.insert("theta", Parameter(1.0));

  _type = "gamma variable";
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
        out, stochbb::gamma(parameter("k").asFloat(), parameter("theta").asFloat(), label().toStdString()));
}

GammaVarNode *
GammaVarNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _type = "compound gamma variable";
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
CompoundGammaVarNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new CompoundGammaVarNode();
}


/* ********************************************************************************************* *
 * Implementation of InvGammaVarNode
 * ********************************************************************************************* */
InvGammaVarNode::InvGammaVarNode(Network *parent)
  : NodeBase(QString("1/")+QChar(0x0393), parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("alpha", Parameter(1.0));
  this->_params.insert("beta", Parameter(1.0));
  _type = "inverse gamma variable";
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
        out, stochbb::invgamma(parameter("alpha").asFloat(), parameter("beta").asFloat(), label().toStdString()));
}

InvGammaVarNode *
InvGammaVarNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _type = "compound inverse gamma variable";
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
CompoundInvGammaVarNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new CompoundInvGammaVarNode();
}


/* ********************************************************************************************* *
 * Implementation of WeibullVarNode
 * ********************************************************************************************* */
WeibullVarNode::WeibullVarNode(Network *parent)
  : NodeBase("Weibull", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("k", Parameter(1.0));
  this->_params.insert("lambda", Parameter(1.0));
  _type = "weibull variable";
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
        out, stochbb::weibull(parameter("k").asFloat(), parameter("lambda").asFloat(), label().toStdString()));
}

WeibullVarNode *
WeibullVarNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _type = "compound weibull variable";
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
CompoundWeibullVarNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new CompoundWeibullVarNode();
}


/* ********************************************************************************************* *
 * Implementation of UniformVarNode
 * ********************************************************************************************* */
UniformVarNode::UniformVarNode(Network *parent)
  : NodeBase("U", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("min", Parameter(0.0));
  this->_params.insert("max", Parameter(1.0));
  _type = "uniform variable";
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
        out, stochbb::uniform(parameter("min").asFloat(), parameter("max").asFloat(), label().toStdString()));
}

UniformVarNode *
UniformVarNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new UniformVarNode();
}


/* ********************************************************************************************* *
 * Implementation of NormalVarNode
 * ********************************************************************************************* */
NormalVarNode::NormalVarNode(Network *parent)
  : NodeBase("N", parent)
{
  this->addSocket(new Socket(QNetSocket::RIGHT, "out", "", this));
  this->_params.insert("mu", Parameter(0.0));
  this->_params.insert("sigma", Parameter(1.0));
  _type = "normal variable";
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
        out, stochbb::normal(parameter("mu").asFloat(), parameter("sigma").asFloat(), label().toStdString()));
}

NormalVarNode *
NormalVarNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _type = "weibull variable";
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
CompoundNormalVarNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _params.insert("graphs", Parameter(int(0)));
  _params.insert("min", Parameter(0.0));
  _params.insert("max", Parameter(1.0));
  _params.insert("steps", Parameter(100));
  _type = "marginal plot";
}

bool
MarginalPlotNode::setParameter(const QString &name, const Parameter &param) {
  if ("graphs" == name) {
    if (param.asInt() < numSockets(QNetSocket::LEFT))
      return false;
    size_t n = numSockets(QNetSocket::LEFT)+1;
    while (param.asInt() > numSockets(QNetSocket::LEFT)) {
      addSocket(new Socket(QNetSocket::LEFT, QString::number(n), QString::number(n), this));
      n++;
    }
  }
  return NodeBase::setParameter(name, param);
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

  size_t nstep = parameter("steps").asInt() > 0 ? parameter("steps").asInt() : 100;
  double tmin = parameter("min").asFloat(), tmax = parameter("max").asFloat();

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
MarginalPlotNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
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
  _params.insert("samples", Parameter(1000));
  _type = "scatter plot";
}

void
ScatterPlotNode::execute(const QHash<Socket *, stochbb::Var> &vartable) {
  stochbb::Var X = vartable[socket("X")];
  stochbb::Var Y = vartable[socket("Y")];

  size_t samples = parameter("samples").asInt() > 0 ? parameter("samples").asInt() : 1000;

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
ScatterPlotNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new ScatterPlotNode();
}


/* ********************************************************************************************* *
 * Implementation of KDEPlotNode
 * ********************************************************************************************* */
KDEPlotNode::KDEPlotNode(Network *parent)
  : OutputNode("KDE Plot", parent)
{
  _params.insert("graphs", Parameter(0));
  _params.insert("samples", Parameter(1000));

  _type = "KDE plot";
}

bool
KDEPlotNode::setParameter(const QString &name, const Parameter &param) {
  if ("graphs" == name) {
    if (param.asInt() < numSockets(QNetSocket::LEFT))
      return false;
    size_t n = numSockets(QNetSocket::LEFT)+1;
    while (param.asInt() > numSockets(QNetSocket::LEFT)) {
      addSocket(new Socket(QNetSocket::LEFT, QString::number(n), QString::number(n), this));
      n++;
    }
  }
  return NodeBase::setParameter(name, param);
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

  size_t nsample = parameter("samples").asInt() > 0 ? parameter("samples").asInt() : 1000;

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
KDEPlotNode::fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable) {
  return new KDEPlotNode();
}


