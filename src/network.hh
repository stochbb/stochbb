#ifndef NETWORK_HH
#define NETWORK_HH

#include "qnetview.hh"
#include <QDomDocument>
#include <QAbstractMessageHandler>
#include "nodes.hh"



class SchemaMessageHandler: public QAbstractMessageHandler
{
  Q_OBJECT

public:
  SchemaMessageHandler(ParserInfo &info, QObject *parent=0);

protected:
  void handleMessage(QtMsgType type, const QString &description,
                     const QUrl &identifier, const QSourceLocation &sourceLocation);

protected:
  ParserInfo &_info;
};


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
  bool load(const QString &file, ParserInfo &info);
  bool load(const QDomDocument &doc, ParserInfo &info);
  bool save();
  bool save(const QString &file);
  QDomDocument serialize() const;

protected:
  QString _filepath;
  void mouseDoubleClickEvent(QMouseEvent *evt);
};

#endif // NETWORK_HH
