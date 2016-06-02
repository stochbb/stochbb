#ifndef NETWORK_HH
#define NETWORK_HH

#include "qnetview.hh"
#include <QDomDocument>

class NodeBase;
class Socket;

class Network : public QNetView
{
  Q_OBJECT

public:
  explicit Network(QWidget *parent = 0);

  virtual bool canConnect(QNetSocket *a, QNetSocket *b);
  virtual void addConnection(QNetSocket *a, QNetSocket *b);
  Socket *findSource(Socket *dest);

  bool hasFilename() const;
  QString filename() const;

  virtual void clear();

signals:
  void editNodeConfig(NodeBase *node);

public slots:
  bool load(const QString &file);
  bool load(const QDomDocument &doc);
  bool save();
  bool save(const QString &file);
  QDomDocument serialize() const;

protected:
  QString _filepath;
  void mouseDoubleClickEvent(QMouseEvent *evt);
};

#endif // NETWORK_HH
