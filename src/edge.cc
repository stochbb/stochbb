#include "edge.hh"
#include "nodes.hh"
#include "network.hh"
#include <QTextStream>


Edge::Edge(Socket *a, Socket *b, Network *net)
  : QNetEdge(a, b, net)
{
  // pass...
}

QDomElement
Edge::serialize(QDomDocument &doc) const {
  QDomElement node = doc.createElement("edge");
  node.setAttribute("srcNode", dynamic_cast<NodeBase *>(_a->parent())->id());
  node.setAttribute("srcSocket", dynamic_cast<Socket *>(_a)->name());
  node.setAttribute("destNode", dynamic_cast<NodeBase *>(_b->parent())->id());
  node.setAttribute("destSocket", dynamic_cast<Socket *>(_b)->name());
  return node;
}

Edge *
Edge::fromXml(const QDomElement &node, const QHash<QString, NodeBase *> &node_table, ParserInfo &info) {
  if (! node.hasAttribute("srcNode")) {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Edge has no source node attribute.";
    info.addError(text);
    return 0;
  }
  if (! node.hasAttribute("srcSocket")) {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Edge has no source socket attribute.";
    info.addError(text);
    return 0;
  }
  if (! node.hasAttribute("destNode")) {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Edge has no destination node attribute.";
    info.addError(text);
    return 0;
  }
  if (! node.hasAttribute("destSocket")) {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Edge has no destination socket attribute.";
    info.addError(text);
    return 0;
  }
  if (! node_table.contains(node.attribute("srcNode"))) {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Unknown source node "
        << node.attribute("srcNode") << ".";
    info.addError(text);
    return 0;
  }
  if (! node_table.contains(node.attribute("destNode"))) {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Unknown destination node "
        << node.attribute("destNode") << ".";
    info.addError(text);
    return 0;
  }

  NodeBase *srcNode = node_table[node.attribute("srcNode")];
  NodeBase *destNode = node_table[node.attribute("destNode")];

  if (! srcNode->hasSocket(node.attribute("srcSocket"))) {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Source node has no socket '"
        << node.attribute("srcSocket") << "'.";
    info.addError(text);
    return 0;
  }
  if (! destNode->hasSocket(node.attribute("destSocket"))) {
    QString text; QTextStream msg(&text);
    msg << "@line " << node.lineNumber() << ": Destination node has no socket '"
        << node.attribute("destSocket") << "'.";
    info.addError(text);
    return 0;
  }

  return new Edge(srcNode->socket(node.attribute("srcSocket")),
                  destNode->socket(node.attribute("destSocket")));
}
