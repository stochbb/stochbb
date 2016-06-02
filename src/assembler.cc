#include "assembler.hh"
#include "network.hh"
#include "nodes.hh"
#include "edge.hh"
#include <QDebug>


/* ********************************************************************************************* *
 * Implementation of Message, MessageBuilder
 * ********************************************************************************************* */
Message::Message(Level level, const QString message)
  : _level(level), _text(message)
{
  // pass...
}

Message::Message(const Message &other)
  : _level(other._level), _text(other._text)
{
  // pass...
}

Message &
Message::operator= (const Message &other) {
  _level = other._level;
  _text  = other._text;
  return *this;
}

Message::Level
Message::level() const {
  return _level;
}

const QString &
Message::text() const {
  return _text;
}

MessageBuilder::MessageBuilder(Message::Level level, Messages &lst)
  : QTextStream(new QString()), _level(level), _lst(lst)
{
  // pass...
}

MessageBuilder::~MessageBuilder() {
  _lst.append(Message(_level, *(this->string())));
  qDebug() << *(this->string());
}


/* ********************************************************************************************* *
 * Implementation of Assembler
 * ********************************************************************************************* */
Assembler::Assembler(Network *net, QHash<Socket *, stochbb::Var> &varTable)
  : _network(net), _queue(), _varTable(varTable), _processedNodes()
{
  // Prepare processing queue
  Network::nodeIterator item = net->nodesBegin();
  for (; item != net->nodesEnd(); item++) {
    if (NodeBase *node = dynamic_cast<NodeBase *>(*item))
      _queue.append(node);
    else
      msgWarn(_messages) << "Cannot cast QNetNode node to NodeBase.";
  }
}

bool
Assembler::assemble() {
  bool progress = false;
  do {
    // Flag if some progress has been made
    // this flag is used to detect unresolvable dependencies
    progress = false;

    // Iterate over queue
    Queue::iterator item = _queue.begin();
    while (item != _queue.end())
    {
      // get node to be processed
      NodeBase *node = *item;

      // if node has been processed already -> remove from queue and continue
      if (_processedNodes.contains(node)) {
        item = _queue.erase(item);
        progress = true;
        continue;
      }

      // Check if node needs pre-processing
      if (needsPreprocessing(node)) {
        if (! preprocess(node))
          return false;
      }

      // Check if all nodes, this node depends on has been processed
      if (processable(node)) {
        // try process node, on error -> fail
        if (! process(node)) {
          msgError(_messages) << "Cannot process node " << node->label() << ".";
          return false;
        }
        // on success, remove node from queue
        item = _queue.erase(item);
        // add node to processed nodes
        _processedNodes.insert(node);
        // flag that some progress has been made
        progress = true;
      } else {
        // if node cannot be processed yet -> continue
        item++;
      }
    }
    // continue as long as we make progress
  } while (progress);

  // If not progress has been made and the queue is empty
  // -> we are stuck, e.g. due to circular dependencies
  if (_queue.size()) {
    msgError(_messages) << "Cannot process network. Some dependencies cannot be resolved.";
    // Stuck...
    return false;
  }

  return true;
}

bool
Assembler::assemble(Network *net, QHash<Socket *, stochbb::Var> &varTable) {
  Assembler ass(net, varTable);
  return ass.assemble();
}

bool
Assembler::assemble(Network *net, QHash<Socket *, stochbb::Var> &varTable, Messages &messages) {
  Assembler ass(net, varTable);
  bool res = ass.assemble();
  messages = ass.messages();
  return res;
}

bool
Assembler::needsPreprocessing(NodeBase *node) {
  if (InhibitionNode *inhib = dynamic_cast<InhibitionNode *>(node)) {
    return (! _varTable.contains(socket(inhib, "Aout"))) ||
        (! _varTable.contains(socket(inhib, "Bout")));
  }
  return false;
}

bool
Assembler::preprocess(NodeBase *node) {
  if (InhibitionNode *inhib = dynamic_cast<InhibitionNode *>(node)) {
    if (! socket(inhib, "Aout"))
      return false;
    if (! _varTable.contains(inhib->socket("Aout")))
      addVariable(inhib->socket("Aout"), stochbb::delta(0));
    if (! _varTable.contains(inhib->socket("Bout")))
      addVariable(inhib->socket("Bout"), stochbb::delta(0));
  }
  return true;
}

bool
Assembler::processable(NodeBase *node) {
  // A node is processable if all its input sockets are connected to processed
  // output sockets (stored in _varTable).
  if (node->numSockets(QNetSocket::LEFT)) {
    // Iterate over all input sockets (on the left side)
    for (size_t i=0; i<node->numSockets(QNetSocket::LEFT); i++) {
      // get socket
      Socket *dest = dynamic_cast<Socket *>(node->socketAt(QNetSocket::LEFT, i));
      if ((! dest) || (! _varTable.contains(dest)))
        return false;
    }
  }
  return true;
}

bool
Assembler::process(NodeBase *node) {
  // Dispatch by node type
  if (StimulusNode *stimulus = dynamic_cast<StimulusNode *>(node)) {
    Socket *out = socket(stimulus, "out");
    if (! out)
      return false;
    addVariable(out, stochbb::delta(stimulus->parameter("time")));
  } else if (DelayNode *delay = dynamic_cast<DelayNode *>(node)) {
    Socket *out = socket(delay, "out");
    stochbb::Var in  = sourceVar(delay, "in");
    if (! out || in.isNull())
      return false;
    addVariable(out, stochbb::delta(delay->parameter("delay")) + in);
  } else if (RandomDelayNode *delay = dynamic_cast<RandomDelayNode *>(node)) {
    Socket *out = socket(delay, "out");
    stochbb::Var in  = sourceVar(delay, "in");
    stochbb::Var d   = sourceVar(delay, "delay");
    if (! out || in.isNull() || d.isNull())
      return false;
    try {
      addVariable(out, d + in);
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing random delay process: "
                          << err.what() << ".";
      return false;
    }
  } else if (GammaProcessNode *gamma = dynamic_cast<GammaProcessNode *>(node)) {
    Socket *out = socket(gamma, "out");
    stochbb::Var in  = sourceVar(gamma, "in");
    if (!out || in.isNull())
      return false;
    try {
      addVariable(out, stochbb::gamma(gamma->parameter("k"), gamma->parameter("theta")) + in);
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing compound Gamma process: "
                          << err.what() << ".";
      return false;
    }
  } else if (CompoundGammaProcessNode *gamma = dynamic_cast<CompoundGammaProcessNode *>(node)) {
    Socket *out = socket(gamma, "out");
    stochbb::Var in  = sourceVar(gamma, "in");
    stochbb::Var k   = sourceVar(gamma, "k");
    stochbb::Var theta = sourceVar(gamma, "theta");
    if (!out || in.isNull() || k.isNull() || theta.isNull())
      return false;
    try {
      addVariable(out, stochbb::gamma(k, theta) + in);
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing compound Gamma process: "
                          << err.what() << ".";
      return false;
    }
  } else if (InvGammaProcessNode *gamma = dynamic_cast<InvGammaProcessNode *>(node)) {
    Socket *out = socket(gamma, "out");
    stochbb::Var in  = sourceVar(gamma, "in");
    if (!out || in.isNull())
      return false;
    try {
      addVariable(out, stochbb::invgamma(gamma->parameter("alpha"), gamma->parameter("beta")) + in);
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing inverse Gamma process: "
                          << err.what() << ".";
      return false;
    }
  } else if (CompoundInvGammaProcessNode *gamma = dynamic_cast<CompoundInvGammaProcessNode *>(node)) {
    Socket *out = socket(gamma, "out");
    stochbb::Var in  = sourceVar(gamma, "in");
    stochbb::Var alpha = sourceVar(gamma, "alpha");
    stochbb::Var beta = sourceVar(gamma, "beta");
    if (!out || in.isNull() || alpha.isNull() || beta.isNull())
      return false;
    try {
      addVariable(out, stochbb::invgamma(alpha, beta) + in);
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing compound inverse Gamma process: "
                          << err.what() << ".";
      return false;
    }
  } else if (MinimumNode *min = dynamic_cast<MinimumNode *>(node)) {
    stochbb::Var X = sourceVar(min, "X");
    stochbb::Var Y = sourceVar(min, "Y");
    Socket *out = socket(min, "out");
    if (!out || X.isNull() || Y.isNull())
      return false;
    try {
      addVariable(out, stochbb::minimum(X, Y));
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing minimum variable: "
                          << err.what() << ".";
      return false;
    }
  } else if (MaximumNode *max = dynamic_cast<MaximumNode *>(node)) {
    stochbb::Var X = sourceVar(max, "X");
    stochbb::Var Y = sourceVar(max, "Y");
    Socket *out = socket(max, "out");
    if (!out || X.isNull() || Y.isNull())
      return false;
    try {
      addVariable(out, stochbb::maximum(X, Y));
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing maximum variable: "
                          << err.what() << ".";
      return false;
    }
  } else if (InhibitionNode *inhb = dynamic_cast<InhibitionNode *>(node)) {
    stochbb::Var X = sourceVar(inhb, "inX");
    stochbb::Var Y = sourceVar(inhb, "inY");
    stochbb::Var A = sourceVar(inhb, "Ain");
    stochbb::Var B = sourceVar(inhb, "Bin");
    Socket *out = socket(inhb, "out");
    if ((!out) || X.isNull() || Y.isNull() || A.isNull() || B.isNull())
      return false;
    try {
      addVariable(out, stochbb::condchain(X, Y, A, B));
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing conditional chain variable: "
                          << err.what() << ".";
      return false;
    }
  } else if (AffineNode *affine = dynamic_cast<AffineNode *>(node)) {
    stochbb::Var in = sourceVar(affine, "in");
    Socket *out = socket(affine, "out");
    if (! out || in.isNull())
      return false;
    addVariable(out, affine->parameter("scale") * in + affine->parameter("shift"));
  } else if (ConstantNode *con = dynamic_cast<ConstantNode *>(node)) {
    Socket *out = socket(con, "out");
    if (! out)
      return false;
    try {
      addVariable(out, stochbb::delta(con->parameter("value")));
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing compound Gamma variable: "
                          << err.what() << ".";
      return false;
    }
  } else if (GammaVarNode *gamma = dynamic_cast<GammaVarNode *>(node)) {
    Socket *out = socket(gamma, "out");
    if (! out)
      return false;
    try {
      addVariable(out, stochbb::gamma(gamma->parameter("k"), gamma->parameter("theta")));
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing Gamma variable: "
                          << err.what() << ".";
      return false;
    }
  } else if (CompoundGammaVarNode *gamma = dynamic_cast<CompoundGammaVarNode *>(node)) {
    stochbb::Var k = sourceVar(gamma, "k");
    stochbb::Var theta = sourceVar(gamma, "theta");
    Socket *out = socket(gamma, "out");
    if (!out || k.isNull() || theta.isNull())
      return false;
    try {
      addVariable(out, stochbb::gamma(k, theta));
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing compound Gamma variable: "
                          << err.what() << ".";
      return false;
    }
  } else if (UniformVarNode *unif = dynamic_cast<UniformVarNode *>(node)) {
    Socket *out = socket(unif, "out");
    if (! out)
      return false;
    try {
      addVariable(out, stochbb::uniform(unif->parameter("min"), unif->parameter("max")));
    } catch (stochbb::Error &err) {
      msgError(_messages) << "Error constructing Uniform variable: "
                          << err.what() << ".";
      return false;
    }
  } else if (dynamic_cast<OutputNode *>(node)) {
    // pass...
  } else {
    msgError(_messages) << "Unknonw type of node " << node->label() << ".";
    return false;
  }

  return true;
}

bool
Assembler::addVariable(Socket *sock, const stochbb::Var &var) {
  if (_varTable.contains(sock)) {
    msgError(_messages) << "Cannot add redefine variable for socket " << sock->name() << ".";
    return false;
  }
  _varTable.insert(sock, var);
  QList<QNetSocket *> dests = _network->findDestinations(sock);
  foreach (QNetSocket *dest, dests) {
    if (Socket *s = dynamic_cast<Socket *>(dest)) {
      _varTable.insert(s, var);
    }
  }
  return true;
}

Socket *
Assembler::socket(NodeBase *node, const QString &name) {
  if (! node->hasSocket(name)) {
    msgError(_messages) << "Node " << node->label() << " has no '" << name << "' socket.";
    return 0;
  }
  return node->socket(name);
}

stochbb::Var
Assembler::sourceVar(NodeBase *node, const QString &name) {
  Socket *sock = node->socket(name);
  if (0 == sock) {
    msgError(_messages) << "No socket '" << name
                        << "' known for node " << node->label() << ".";
    return stochbb::Var();
  }

  if (! _varTable.contains(sock)) {
    msgError(_messages) << "Source connected to '" << sock << "' socket of node "
                        << node->label() << " not processed yet!";
    return stochbb::Var();
  }

  return _varTable[sock];
}

const Messages &
Assembler::messages() const {
  return _messages;
}
