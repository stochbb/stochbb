#include "network.hh"
#include "nodes.hh"
#include "edge.hh"

#include <QDomNodeList>
#include <QDomDocument>
#include <QFile>
#include <QFileInfo>


Network::Network(QWidget *parent)
  : QNetView(parent)
{
  // pass...
}

bool
Network::canConnect(QNetSocket *a, QNetSocket *b) {
  if (! QNetView::canConnect(a, b))
    return false;

  // Connect only OUT -> IN
  if ((QNetSocket::RIGHT != a->side()) || (QNetSocket::LEFT != b->side()))
    return false;

  // Check if input is not used yet
  foreach (QNetEdge *edge, _edges) {
    if (edge->dest() == b)
      return false;
  }

  return dynamic_cast<Socket *>(a) && dynamic_cast<Socket *>(b);
}

void
Network::addConnection(QNetSocket *a, QNetSocket *b) {
  addEdge(new Edge(dynamic_cast<Socket *>(a), dynamic_cast<Socket *>(b)));
}

Socket *
Network::findSource(Socket *dest) {
  QList<QNetSocket *> srcs = this->findSources(dest);
  if (1 != srcs.size())
    return 0;
  return dynamic_cast<Socket *>(srcs.front());
}

void
Network::mouseDoubleClickEvent(QMouseEvent *evt) {
  QNetView::mouseDoubleClickEvent(evt);

  if (_selectedNode && dynamic_cast<NodeBase *>(_selectedNode))
    emit editNodeConfig(dynamic_cast<NodeBase *>(_selectedNode));
}

bool
Network::hasFilename() const {
  return ! _filepath.isEmpty();
}

QString
Network::filename() const {
  QFileInfo info(_filepath);
  return info.fileName();
}

void
Network::clear() {
  QNetView::clear();
  _filepath.clear();
}

bool
Network::load(const QString &file) {
  QFile fd(file);
  if (! fd.open(QIODevice::ReadOnly))
    return false;

  QDomDocument doc;
  if (! doc.setContent(& fd))
    return false;

  bool success = this->load(doc);

  if (success)
    _filepath = file;

  return success;
}

bool
Network::save() {
  if (_filepath.isEmpty())
    return false;
  return save(_filepath);
}

bool
Network::save(const QString &file) {
  QDomDocument doc = this->serialize();

  QFile fd(file);
  if (! fd.open(QIODevice::WriteOnly))
    return false;

  fd.write(doc.toString().toUtf8());
  fd.flush();
  fd.close();

  _filepath = file;
  return true;
}


bool
Network::load(const QDomDocument &doc) {
  clear();

  QDomElement root = doc.documentElement();
  QHash<QString, NodeBase *> node_table;

  // Iterate over all "node" children
  QDomElement node = root.firstChildElement("node");
  for (; ! node.isNull(); node = node.nextSiblingElement("node")) {
    NodeBase *obj = NodeBase::fromXml(node);
    if (! obj)
      return false;
    node_table.insert(node.attribute("id"), obj);
    this->addNode(obj);
  }

  // Iterate over all "edge" children
  QDomElement edge = root.firstChildElement("edge");
  for (; ! edge.isNull(); edge = edge.nextSiblingElement("edge")) {
    Edge *obj = Edge::fromXml(edge, node_table);
    if (! obj)
      return false;
    this->addEdge(obj);
  }

  return true;
}

QDomDocument
Network::serialize() const {
  QDomDocument doc;
  QDomElement root = doc.createElement("net");
  foreach (QNetNode *obj, _nodes) {
    if (NodeBase *node = dynamic_cast<NodeBase *>(obj)) {
      root.appendChild(node->serialize(doc));
    }
  }
  foreach (QNetEdge *obj, _edges) {
    if (Edge *edge = dynamic_cast<Edge *>(obj)) {
      root.appendChild(edge->serialize(doc));
    }
  }
  doc.appendChild(root);
  return doc;
}
