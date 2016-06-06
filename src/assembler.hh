#ifndef ASSEMBLER_HH
#define ASSEMBLER_HH

#include <stochbb/api.hh>
#include <QHash>
#include <QList>
#include <QSet>
#include <QTextStream>

#include "nodes.hh"

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
  bool hasVariable(Socket *sock) const;
  bool addVariable(Socket *sock, const stochbb::Var &var);
  Socket *socket(const NodeBase *node, const QString &name);
  stochbb::Var sourceVar(const NodeBase *node, const QString &name);

public:
  static bool assemble(Network *net, QHash<Socket *, stochbb::Var> &varTable);
  static bool assemble(Network *net, QHash<Socket *, stochbb::Var> &varTable, Messages &messages);

protected:
  typedef QList<NodeBase *> Queue;

protected:
  Assembler(Network *net, QHash<Socket *, stochbb::Var> &varTable);
  bool assemble();

  const Messages &messages() const;

protected:
  Network *_network;
  Queue    _queue;
  QHash<Socket *, stochbb::Var> &_varTable;
  QSet<NodeBase *> _processedNodes;
  Messages _messages;
};

#endif // ASSEMBLER_HH
