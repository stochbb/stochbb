#ifndef EDGE_HH
#define EDGE_HH

#include "qnetview.hh"
#include <QDomElement>

class Socket;
class Network;
class NodeBase;
class ParserInfo;

class Edge : public QNetEdge
{
  Q_OBJECT

public:
  explicit Edge(Socket *a, Socket *b, Network *net=0);

  QDomElement serialize(QDomDocument &doc) const;

public:
  static Edge *fromXml(const QDomElement &node, const QHash<QString, NodeBase *> &node_table, ParserInfo &info);
};

#endif // EDGE_HH
