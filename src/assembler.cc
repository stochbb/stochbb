#include "assembler.hh"
#include "network.hh"
#include "nodes.hh"
#include "edge.hh"


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
      if (node->needsPreprocessing(*this)) {
        if (! node->preprocess(*this))
          return false;
      }

      // Check if all nodes, this node depends on has been processed
      if (node->processable(*this)) {
        // try process node, on error -> fail
        try {
          if(! node->assemble(*this)) {
            msgError(_messages) << "Cannot process node " << node->label() << ".";
            return false;
          }
        } catch (stochbb::Error &error) {
          msgError(_messages) << "Cannot process node " << node->label() << ": "
                              << error.str().c_str() << ".";
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
    msgError(_messages) << "Cannot process network. Some dependencies cannot be resolved. "
                        << "This is usually due to unconnected inputs or circular dependencies "
                        << "between nodes.";
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
Assembler::hasVariable(Socket *sock) const {
  return _varTable.contains(sock);
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
Assembler::socket(const NodeBase *node, const QString &name) {
  if (! node->hasSocket(name)) {
    msgError(_messages) << "Node " << node->label() << " has no '" << name << "' socket.";
    return 0;
  }
  return node->socket(name);
}

stochbb::Var
Assembler::sourceVar(const NodeBase *node, const QString &name) {
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
