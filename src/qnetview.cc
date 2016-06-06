#include "qnetview.hh"
#include <QPainter>
#include <QFontMetrics>
#include <QPaintEvent>
#include <QDebug>
#include <QTransform>

#define SOCKET_SIZE 7
#define SOCKET_LABEL_PADDING 3
#define EDGE_CURVE_DIST 20


/* ********************************************************************************************* *
 * Implementation of QNetNode
 * ********************************************************************************************* */
QNetNode::QNetNode(const QString &label, QNetView *parent)
  : QObject(parent), _borderPen(Qt::gray, 3), _altBorderPen(Qt::blue, 3),
    _backgroundBrush(Qt::white), _margin(5), _socketMargin(3),
    _selected(false), _position(0,0), _labelFont("sans", 12, QFont::Bold), _label(label),
    _description()
{
  updateLayout();
}

const QString &
QNetNode::label() const {
  return _label;
}

void
QNetNode::setLabel(const QString &label) {
  _label = label;
  updateLayout();
}

bool
QNetNode::hasDescription() const {
  return ! _description.isEmpty();
}

const QString &
QNetNode::description() const {
  return _description;
}

void
QNetNode::setDescription(const QString &text) {
  _description = text;
}

const QPoint &
QNetNode::position() const {
  return _position;
}

const QSize &
QNetNode::size() const {
  return _size;
}

void
QNetNode::setPosition(const QPoint &pos) {
  _position = pos;
  updateLayout();
}

void
QNetNode::setPosition(int x, int y) {
  setPosition(QPoint(x,y));
}

size_t
QNetNode::numSockets(QNetSocket::Side side) const {
  switch (side) {
    case QNetSocket::LEFT: return _leftSockets.size();
    case QNetSocket::RIGHT: return _rightSockets.size();
  }
  return 0;
}

void
QNetNode::addSocket(QNetSocket *socket) {
  switch (socket->side()) {
    case QNetSocket::LEFT:
      _leftSockets.append(socket);
      break;
    case QNetSocket::RIGHT:
      _rightSockets.append(socket);
      break;
  }
  socket->setParent(this);
  updateLayout();
}

QNetSocket *
QNetNode::socketAt(const QPoint &pos) const {
  if ((pos.x() <= _position.x()) && (pos.x() >= (_position.x()-SOCKET_SIZE))) {
    foreach (QNetSocket *socket, _leftSockets) {
      QRect rect(socket->position(), socket->size());
      if (rect.contains(pos)) {
        return socket;
      }
    }
  } else if ((pos.x() >= (_position.x()+_size.width())) && (pos.x() <= (_position.x()+_size.width()+SOCKET_SIZE))) {
    foreach (QNetSocket *socket, _rightSockets) {
      QRect rect(socket->position(), socket->size());
      if (rect.contains(pos)) {
        return socket;
      }
    }
  }
  return 0;
}

QNetSocket *
QNetNode::socketAt(QNetSocket::Side side, size_t i) const {
  if (QNetSocket::LEFT == side) {
    if (i >= _leftSockets.size())
      return 0;
    return _leftSockets.at(i);
  }
  if (i >= _rightSockets.size())
    return 0;
  return _rightSockets.at(i);
}

bool
QNetNode::selected() const {
  return _selected;
}

void
QNetNode::select(bool selected) {
  _selected = selected;
}

QRect
QNetNode::boundingRect() const {
  QPoint pos = position();
  QSize  size = this->size();
  pos -= QPoint(_borderPen.width()/2,_borderPen.width()/2);
  size += QSize(_borderPen.width(),_borderPen.width());
  QRect rect(pos, size);
  foreach(QNetSocket *socket, _leftSockets) {
    rect = rect.united(socket->boundingRect());
  }
  foreach(QNetSocket *socket, _rightSockets) {
    rect = rect.united(socket->boundingRect());
  }
  return rect;
}

void
QNetNode::updateLayout() {
  int ls_size = 2*_margin, rs_size = 2*_margin;
  if (_leftSockets.size())
    ls_size += _leftSockets.size()*SOCKET_SIZE + (_leftSockets.size()-1)*_socketMargin;
  if (_rightSockets.size())
    rs_size += _rightSockets.size()*SOCKET_SIZE + (_rightSockets.size()-1)*_socketMargin;

  int label_top_pad = _margin, label_bottom_pad = _margin, label_left_pad = _margin,
      label_right_pad = _margin;
  int max_left_sock_width = _margin, max_right_sock_width = _margin;
  // update label padding on socket labels
  foreach (QNetSocket *socket, _leftSockets) {
    max_left_sock_width = std::max(max_left_sock_width, socket->size().width());
  }
  foreach (QNetSocket *socket, _rightSockets) {
    max_right_sock_width = std::max(max_right_sock_width, socket->size().width());
  }
  label_left_pad += max_left_sock_width;
  label_right_pad += max_right_sock_width;

  QFontMetrics fm(_labelFont);
  int label_height = fm.height() + label_top_pad + label_bottom_pad;
  int label_width  = fm.width(_label) + label_left_pad + label_right_pad;
  int sock_height  = std::max(ls_size, rs_size);
  int height = std::max(label_height, sock_height);
  int width  = label_width;
  _size = QSize(width, height);

  // compute socket positions
  _leftSocketOffset = _margin;
  if (height > ls_size)
    _leftSocketOffset += (height-ls_size)/2;
  _rightSocketOffset = _margin;
  if (height > rs_size)
    _rightSocketOffset += (height-rs_size)/2;

  // compute label position
  _labelPos = QPoint(_position.x()+label_left_pad, _position.y()+label_top_pad + fm.ascent());
  if (height > label_height)
    _labelPos += QPoint(0, (height - label_height)/2);

  // Update positions
  int offset = _leftSocketOffset;
  foreach (QNetSocket *socket, _leftSockets) {
    socket->setPosition(_position.x()-SOCKET_SIZE, _position.y()+offset);
    offset += (SOCKET_SIZE + _socketMargin);
  }
  offset = _rightSocketOffset;
  foreach (QNetSocket *socket, _rightSockets) {
    socket->setPosition(_position.x()+width, _position.y()+offset);
    offset += (SOCKET_SIZE + _socketMargin);
  }
}

void
QNetNode::paint(QPainter &painter) {
  if (_selected)
    painter.setPen(_altBorderPen);
  else
    painter.setPen(_borderPen);
  painter.setBrush(_backgroundBrush);
  painter.drawRect(_position.x(), _position.y(), _size.width(), _size.height());

  // draw sockets
  foreach (QNetSocket *socket, _leftSockets) {
    socket->paint(painter);
  }
  foreach (QNetSocket *socket, _rightSockets) {
    socket->paint(painter);
  }

  // draw label
  painter.setFont(_labelFont);
  painter.setPen(Qt::black);
  painter.drawText(_labelPos, _label);
}


/* ********************************************************************************************* *
 * Implementation of QNetSocket
 * ********************************************************************************************* */
QNetSocket::QNetSocket(Side loc, QNetNode *node)
  : QObject(node), _side(loc), _label(), _borderPen(Qt::gray, 2), _labelFont("sans", 8),
    _position(0,0), _size(SOCKET_SIZE, SOCKET_SIZE)
{
  // pass...
}

QNetSocket::QNetSocket(Side loc, const QString &label, QNetNode *node)
  : QObject(node), _side(loc), _label(label), _borderPen(Qt::gray, 2), _labelFont("sans", 8),
    _position(0,0), _size(SOCKET_SIZE, SOCKET_SIZE)
{
  QFontMetrics fm(_labelFont);
  _size += QSize(fm.width(_label)+SOCKET_LABEL_PADDING, 0);
}

const QPoint &
QNetSocket::position() const {
  return _position;
}

const QNetSocket::Side &
QNetSocket::side() const {
  return _side;
}

const QPoint &
QNetSocket::anchor() const {
  return _anchor;
}

const QSize &
QNetSocket::size() const {
  return _size;
}

bool
QNetSocket::hasLabel() const {
  return ! _label.isEmpty();
}

const QString &
QNetSocket::label() const {
  return _label;
}

void
QNetSocket::setPosition(const QPoint &position) {
  _position = position;
  switch (_side) {
    case LEFT:
      _anchor = QPoint(_position.x()-_borderPen.width()/2, _position.y()+SOCKET_SIZE/2);
      break;
    case RIGHT:
      _anchor = QPoint(_position.x()+SOCKET_SIZE+_borderPen.width()/2, _position.y()+SOCKET_SIZE/2);
      break;
  }
}

void
QNetSocket::setPosition(int x, int y) {
  setPosition(QPoint(x,y));
}

QRect
QNetSocket::boundingRect() const {
  QPoint pos = position();
  QSize size = this->size();
  pos -= QPoint(_borderPen.width()/2, _borderPen.width()/2);
  size += QSize(_borderPen.width(), _borderPen.width());
  return QRect(pos, size);
}

void
QNetSocket::paint(QPainter &painter) {
  QPen pen(Qt::gray);
  pen.setWidth(2);
  painter.setPen(pen);
  painter.drawRect(_position.x(), _position.y(), SOCKET_SIZE, SOCKET_SIZE);

  if (! _label.isEmpty()) {
    QFontMetrics fm(_labelFont);
    int x = _position.x();
    int width = fm.width(_label);
    if (LEFT == _side)
      x += SOCKET_SIZE+SOCKET_LABEL_PADDING;
    else
      x -= SOCKET_LABEL_PADDING+width;
    painter.setFont(_labelFont);
    painter.drawText(x, _position.y()+fm.ascent()-fm.strikeOutPos(), _label);
  }
}

/* ********************************************************************************************* *
 * Implementation of QNetEdge
 * ********************************************************************************************* */
QNetEdge::QNetEdge(QNetSocket *a, QNetSocket *b, QNetView *parent)
  : QObject(parent), _pen(Qt::darkBlue, 3), _altPen(Qt::blue, 3), _a(a), _b(b), _selected(false)
{
  connect(a, SIGNAL(destroyed(QObject*)), this, SLOT(deleteLater()));
  connect(b, SIGNAL(destroyed(QObject*)), this, SLOT(deleteLater()));
}

bool
QNetEdge::selected() const {
  return _selected;
}

bool
QNetEdge::contains(const QPoint &pos) const {
  QRect fuzzy(pos.x()-3,pos.y()-3, 6,6);
  return this->path().intersects(fuzzy);
}

QPainterPath
QNetEdge::path() const {
  QPoint c1, c2;
  switch (_a->side()) {
    case QNetSocket::LEFT:
      c1 = _a->anchor() - QPoint(EDGE_CURVE_DIST, 0);
      break;
    case QNetSocket::RIGHT:
      c1 = _a->anchor() + QPoint(EDGE_CURVE_DIST, 0);
      break;
  }
  switch (_b->side()) {
    case QNetSocket::LEFT:
      c2 = _b->anchor() - QPoint(EDGE_CURVE_DIST, 0);
      break;
    case QNetSocket::RIGHT:
      c2 = _b->anchor() + QPoint(EDGE_CURVE_DIST, 0);
      break;
  }

  QPainterPath path(_a->anchor());
  path.cubicTo(c1, c2, _b->anchor());
  return path;
}

bool
QNetEdge::isConnectedWith(QNetNode *node) const {
  return (_a->parent() == node) ||
      (_b->parent() == node);
}

QNetSocket *
QNetEdge::src() const {
  return _a;
}

QNetSocket *
QNetEdge::dest() const {
  return _b;
}

QRect
QNetEdge::boundingRect() const {
  QPainterPath path = this->path();
  return path.boundingRect().toRect();
}

void
QNetEdge::paint(QPainter &painter) {
  QPainterPath path = this->path();

  if (_selected)
    painter.setPen(_altPen);
  else
    painter.setPen(_pen);

  painter.setBrush(Qt::NoBrush);
  painter.drawPath(path);
}

void
QNetEdge::select(bool selected) {
  _selected = selected;
}


/* ********************************************************************************************* *
 * Implementation of QNetView
 * ********************************************************************************************* */
QNetView::QNetView(QWidget *parent)
  : QWidget(parent), _scale(1), _modified(false), _dragging(0), _dragPos(0,0), _connecting(0),
    _selectedNode(0), _selectedEdge(0)
{
  updateLayout();
}

double
QNetView::scale() const {
  return _scale;
}

void
QNetView::setScale(double scale) {
  _scale = scale;
  updateLayout();
}

void
QNetView::addNode(QNetNode *node) {
  node->setParent(this);
  _nodes.append(node);
  connect(node, SIGNAL(destroyed(QObject*)), this, SLOT(itemDestroyed(QObject*)));
  setModified();
  updateLayout();
}

void
QNetView::remNode(QNetNode *node) {
  if (! _nodes.contains(node))
    return;
  _nodes.removeAll(node);
  node->deleteLater();
}

void
QNetView::clear() {
  bool mod = _nodes.size();
  while (_nodes.size()) {
    this->remNode(_nodes.first());
  }
  setModified(mod);
}

QNetView::nodeIterator
QNetView::nodesBegin() {
  return _nodes.begin();
}

QNetView::nodeIterator
QNetView::nodesEnd() {
  return _nodes.end();
}

bool
QNetView::canConnect(QNetSocket *a, QNetSocket *b) {
  return a != b;
}

void
QNetView::addConnection(QNetSocket *a, QNetSocket *b) {
  this->addEdge(new QNetEdge(a, b, this));
}

void
QNetView::addEdge(QNetEdge *edge) {
  edge->setParent(this);
  _edges.append(edge);
  connect(edge, SIGNAL(destroyed(QObject*)), this, SLOT(itemDestroyed(QObject*)));
  setModified();
  updateLayout();
}

void
QNetView::remEdge(QNetEdge *edge) {
  if (! _edges.contains(edge))
    return;
  _edges.removeAll(edge);
  edge->deleteLater();
}

QNetView::edgeIterator
QNetView::edgesBegin() {
  return _edges.begin();
}

QNetView::edgeIterator
QNetView::edgesEnd() {
  return _edges.end();
}

bool
QNetView::selected() const {
  return (_selectedNode || _selectedEdge);
}

QNetNode *
QNetView::selectedNode() const {
  return _selectedNode;
}

QNetEdge *
QNetView::selectedEdge() const {
  return _selectedEdge;
}

QList<QNetSocket *>
QNetView::findSources(QNetSocket *dest) {
  QList<QNetSocket *> sockets;
  foreach (QNetEdge *edge, _edges) {
    if (edge->dest() == dest)
      sockets.append(edge->src());
  }
  return sockets;
}

QList<QNetSocket *>
QNetView::findDestinations(QNetSocket *src) {
  QList<QNetSocket *> sockets;
  foreach (QNetEdge *edge, _edges) {
    if (edge->src() == src)
      sockets.append(edge->dest());
  }
  return sockets;
}

bool
QNetView::isModified() const{
  return _modified;
}

QRect
QNetView::boundingRect() const {
  if (0 == _nodes.size())
    return QRect();

  QRect rect;
  foreach(QNetNode *node, _nodes) {
    rect = rect.united(node->boundingRect());
  }
  foreach(QNetEdge *edge, _edges) {
    rect = rect.united(edge->boundingRect());
  }

  return rect;
}

void
QNetView::setModified(bool mod) {
  bool ismod = _modified ^ mod;
  _modified = mod;
  if (ismod)
    emit modified();
}

void
QNetView::itemDestroyed(QObject *item) {
  if (_nodes.contains((QNetNode *) item)) {
    _nodes.removeAll((QNetNode *) item);
    if (_selectedNode == (QNetNode *) item)
      _selectedNode = 0;
  } else if (_edges.contains((QNetEdge *) item)) {
    _edges.removeAll((QNetEdge *) item);
    if (_selectedEdge == (QNetEdge *) item)
      _selectedEdge = 0;
  } else {
    return;
  }
  setModified();
  updateLayout();
}


void
QNetView::updateLayout() {
  QRect bb(QPoint(), this->minimumSize());
  bb = bb.united(QRect(QPoint(0,0), this->size()));

  foreach (QNetNode *node, _nodes) {
    bb = bb.united(
          QRect(node->position()-QPoint(SOCKET_SIZE,0), node->size()+QSize(SOCKET_SIZE,0)));
  }
  foreach (QNetEdge *edge, _edges) {
    bb = bb.united(edge->path().boundingRect().toRect());
  }

  this->setMinimumSize(bb.right()*_scale,bb.bottom()*_scale);

  update();
}

void
QNetView::paint(QPainter &painter) {
  foreach (QNetNode *node, _nodes) {
    node->paint(painter);
  }
  foreach (QNetEdge *edge, _edges) {
    edge->paint(painter);
  }

  if (_connecting) {
    QPainterPath path(_connecting->anchor());
    QPoint c1,c2;
    switch (_connecting->side()) {
      case QNetSocket::LEFT:
        c1 = _connecting->anchor() - QPoint(30, 0);
        c2 = _dragPos + QPoint(30, 0);
        break;
      case QNetSocket::RIGHT:
        c1 = _connecting->anchor() + QPoint(30, 0);
        c2 = _dragPos - QPoint(30, 0);
        break;
    }
    path.cubicTo(c1, c2, _dragPos);
    QPen pen(Qt::black, 2, Qt::DotLine);
    painter.setPen(pen);
    painter.setBrush(Qt::NoBrush);
    painter.drawPath(path);
  }
}

void
QNetView::paintEvent(QPaintEvent *evt) {
  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);
  painter.fillRect(evt->rect(), Qt::white);

  QTransform scale; scale.scale(_scale, _scale);
  painter.setTransform(scale);

  paint(painter);
}

void
QNetView::mousePressEvent(QMouseEvent *evt) {
  QWidget::mousePressEvent(evt);

  if (_selectedNode) {
    _selectedNode->select(false);
    _selectedNode=0;
  }
  if (_selectedEdge) {
    _selectedEdge->select(false);
    _selectedEdge=0;
  }

  if (Qt::LeftButton == evt->button()) {
    QList<QNetNode *>::reverse_iterator node = _nodes.rbegin();
    for (; node != _nodes.rend(); node++) {
      QRect bb((*node)->position(), (*node)->size());
      if (bb.contains(evt->pos()/_scale)) {
        _dragging = _selectedNode = *node;
        _selectedNode->select(true);
        _dragPos = (*node)->position()-evt->pos()/_scale;
        break;
      } else if (QNetSocket *socket = (*node)->socketAt(evt->pos()/_scale)) {
        _connecting = socket;
        _dragPos = evt->pos()/_scale;
        break;
      }
    }
    if (! _selectedNode) {
      foreach (QNetEdge *edge, _edges) {
        if (edge->contains(evt->pos()/_scale)) {
          _selectedEdge = edge;
          _selectedEdge->select(true);
          break;
        }
      }
    }
  } else if (Qt::RightButton == evt->button()) {
    foreach (QNetNode *node, _nodes) {
      QRect bb(node->position(), node->size());
      if (bb.contains(evt->pos()/_scale)) {
        _selectedNode = node; _selectedNode->select(true);
        break;
      }
    }
    if (! _selectedNode) {
      foreach (QNetEdge *edge, _edges) {
        if (edge->contains(evt->pos()/_scale)) {
          _selectedEdge = edge;
          _selectedEdge->select(true);
          break;
        }
      }
    }
  }

  update();
}

void
QNetView::mouseMoveEvent(QMouseEvent *evt) {
  QWidget::mouseMoveEvent(evt);

  if ((0 >= evt->pos().x()) || (0 >= evt->pos().y()))
    return;

  if (_dragging) {
    _dragging->setPosition(evt->pos()/_scale+_dragPos);
    setModified();
    updateLayout();
    update();
  }
  if (_connecting) {
    _dragPos = evt->pos()/_scale;
    update();
  }
}

void
QNetView::mouseReleaseEvent(QMouseEvent *evt) {
  QWidget::mouseReleaseEvent(evt);
  _dragging = 0;

  if (_connecting) {
    foreach (QNetNode *node, _nodes) {
      if (QNetSocket *socket = node->socketAt(evt->pos()/_scale)) {
        if (this->canConnect(_connecting, socket))
          this->addConnection(_connecting, socket);
        break;
      }
    }
    _connecting = 0;
    update();
  }
}
