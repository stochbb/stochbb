#ifndef QNETVIEW_HH
#define QNETVIEW_HH

#include <QWidget>
#include <QSet>
#include <QList>
#include <QFont>
#include <QPen>
#include <QBrush>


// Forward declarations
class QNetView;
class QNetNode;
class QNetEdge;
class QNetSocket;


class QNetSocket: public QObject
{
  Q_OBJECT

public:
  typedef enum {
    LEFT, RIGHT
  } Side;

public:
  explicit QNetSocket(Side loc, QNetNode *node);
  explicit QNetSocket(Side loc, const QString &label, QNetNode *node);

  const QPoint &position() const;
  const Side &side() const;
  const QSize &size() const;
  const QPoint &anchor() const;

  void setPosition(const QPoint &position);
  void setPosition(int x, int y);

  bool hasLabel() const;
  const QString &label() const;

  QRect boundingRect() const;

public slots:
  void paint(QPainter &painter);

protected:
  Side     _side;
  QString  _label;
  QPen     _borderPen;
  QFont    _labelFont;

  QPoint   _position;
  QSize    _size;
  QPoint   _anchor;
};


class QNetNode: public QObject
{
  Q_OBJECT

public:
  explicit QNetNode(const QString &label="", QNetView *parent=0);

  const QString &label() const;
  void setLabel(const QString &label);

  bool hasDescription() const;
  const QString &description() const;
  void setDescription(const QString &text);

  const QPoint &position() const;
  const QSize &size() const;

  void setPosition(const QPoint &pos);
  void setPosition(int x, int y);

  size_t numSockets(QNetSocket::Side side) const;
  virtual void addSocket(QNetSocket *socket);
  QNetSocket *socketAt(const QPoint &pos) const;
  QNetSocket *socketAt(QNetSocket::Side side, size_t i) const;

  bool selected() const;

  QRect boundingRect() const;

public slots:
  void paint(QPainter &painter);

protected slots:
  void updateLayout();
  void select(bool selected);

protected:
  QPen    _borderPen;
  QPen    _altBorderPen;
  QBrush  _backgroundBrush;
  int     _margin;
  int     _socketMargin;

  bool    _selected;
  QPoint  _position;
  QSize   _size;

  QPoint  _labelPos;
  QFont   _labelFont;
  QString _label;

  QString _description;

  int     _leftSocketOffset;
  int     _rightSocketOffset;
  QList<QNetSocket *> _leftSockets;
  QList<QNetSocket *> _rightSockets;

  friend class QNetView;
};


class QNetEdge: public QObject
{
  Q_OBJECT

public:
  explicit QNetEdge(QNetSocket *a, QNetSocket *b, QNetView *parent=0);

  bool isConnectedWith(QNetNode *node) const;
  QNetSocket *src() const;
  QNetSocket *dest() const;

  bool selected() const;
  bool contains(const QPoint &pos) const;
  QPainterPath path() const;

  QRect boundingRect() const;

public slots:
  void paint(QPainter &painter);

protected slots:
  void select(bool selected);

protected:
  QPen _pen;
  QPen _altPen;
  QNetSocket *_a;
  QNetSocket *_b;
  bool _selected;

  friend class QNetView;
};


class QNetView: public QWidget
{
  Q_OBJECT

public:
  typedef QList<QNetNode *>::iterator nodeIterator;
  typedef QList<QNetEdge *>::iterator edgeIterator;

public:
  explicit QNetView(QWidget *parent = 0);

  double scale() const;
  void setScale(double scale);

  virtual void addNode(QNetNode *node);
  virtual void remNode(QNetNode *node);
  virtual void clear();
  nodeIterator nodesBegin();
  nodeIterator nodesEnd();

  virtual bool canConnect(QNetSocket *a, QNetSocket *b);
  virtual void addConnection(QNetSocket *a, QNetSocket *b);
  virtual void addEdge(QNetEdge *edge);
  virtual void remEdge(QNetEdge *edge);
  edgeIterator edgesBegin();
  edgeIterator edgesEnd();

  bool selected() const;
  QNetNode *selectedNode() const;
  QNetEdge *selectedEdge() const;

  QList<QNetSocket *> findSources(QNetSocket *dest);
  QList<QNetSocket *> findDestinations(QNetSocket *src);

  bool isModified() const;

  QRect boundingRect() const;

public slots:
  void setModified(bool modified=true);
  void paint(QPainter &painter);

signals:
  void modified();

protected slots:
  void itemDestroyed(QObject *item);

protected:
  void updateLayout();

  void paintEvent(QPaintEvent *evt);
  void mousePressEvent(QMouseEvent *evt);
  void mouseMoveEvent(QMouseEvent *evt);
  void mouseReleaseEvent(QMouseEvent *evt);

protected:
  double _scale;
  QList<QNetNode*> _nodes;
  QList<QNetEdge*> _edges;
  bool _modified;

  QNetNode *_dragging;
  QPoint   _dragPos;
  QNetSocket *_connecting;
  QNetNode *_selectedNode;
  QNetEdge *_selectedEdge;
};

#endif // QNETVIEW_HH
