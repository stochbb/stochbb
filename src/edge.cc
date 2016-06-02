#include "edge.hh"
#include "nodes.hh"
#include "network.hh"

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
Edge::fromXml(const QDomElement &node, const QHash<QString, NodeBase *> &node_table) {
  if (! node.hasAttribute("srcNode"))
    return 0;
  if (! node.hasAttribute("srcSocket"))
    return 0;
  if (! node.hasAttribute("destNode"))
    return 0;
  if (! node.hasAttribute("destSocket"))
    return 0;
  if (! node_table.contains(node.attribute("srcNode")))
    return 0;
  if (! node_table.contains(node.attribute("destNode")))
    return 0;

  NodeBase *srcNode = node_table[node.attribute("srcNode")];
  NodeBase *destNode = node_table[node.attribute("destNode")];
  if (! srcNode->hasSocket(node.attribute("srcSocket")))
    return 0;
  if (! destNode->hasSocket(node.attribute("destSocket")))
    return 0;

  return new Edge(srcNode->socket(node.attribute("srcSocket")),
                  destNode->socket(node.attribute("destSocket")));
}
