#ifndef ASSEMBLER_HH
#define ASSEMBLER_HH

#include "lib/api.hh"
#include <QHash>
#include <QList>
#include <QSet>
#include <QTextStream>


class Socket;
class NodeBase;
class Network;

class Message
{
public:
  typedef enum {
    INFO, WARNING, CRITICAL
  } Level;

public:
  Message(Level level, const QString message);
  Message(const Message &other);
  Message &operator=(const Message &other);

  Level level() const;
  const QString &text() const;

protected:
  Level _level;
  QString _text;
};

typedef QList<Message> Messages;

class MessageBuilder: public QTextStream
{
public:
  MessageBuilder(Message::Level level, Messages &lst);
  virtual ~MessageBuilder();

protected:
  Message::Level _level;
  Messages &_lst;
};

#define msgInfo(lst) MessageBuilder(Message::INFO, lst)
#define msgWarn(lst) MessageBuilder(Message::WARNING, lst)
#define msgError(lst) MessageBuilder(Message::CRITICAL, lst)


class Assembler
{
public:
  static bool assemble(Network *net, QHash<Socket *, stochbb::Var> &varTable);
  static bool assemble(Network *net, QHash<Socket *, stochbb::Var> &varTable, Messages &messages);
protected:
  typedef QList<NodeBase *> Queue;

protected:
  Assembler(Network *net, QHash<Socket *, stochbb::Var> &varTable);
  bool assemble();
  bool needsPreprocessing(NodeBase *node);
  bool preprocess(NodeBase *node);
  bool processable(NodeBase *node);
  bool process(NodeBase *node);
  bool addVariable(Socket *sock, const stochbb::Var &var);
  Socket *socket(NodeBase *node, const QString &name);
  stochbb::Var sourceVar(NodeBase *node, const QString &name);

  const Messages &messages() const;

protected:
  Network *_network;
  Queue    _queue;
  QHash<Socket *, stochbb::Var> &_varTable;
  QSet<NodeBase *> _processedNodes;
  Messages _messages;
};

#endif // ASSEMBLER_HH
