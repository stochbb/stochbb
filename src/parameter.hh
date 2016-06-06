#ifndef PARAMETER_HH
#define PARAMETER_HH

#include <QString>
#include <QVariant>
#include <QWidget>
#include <QDomDocument>
#include <QDomElement>
#include <QDialog>


class QLineEdit;
class QTextEdit;
class NodeBase;


class ParserInfo
{
public:
  typedef enum {
    OK = 0,
    WARNING,
    ERROR
  } State;

public:
  ParserInfo();

  State state() const;
  QStringList messages() const;
  void addMessage(State state, const QString &msg);
  void addInfo(const QString &msg);
  void addWarning(const QString &msg);
  void addError(const QString &msg);

protected:
  State _state;
  QList< QPair<State, QString> > _messages;
};


class Parameter
{
public:
  typedef enum {
    BOOL, INTEGER, FLOAT, STRING
  } Type;

public:
  Parameter();
  explicit Parameter(bool val);
  explicit Parameter(int val);
  explicit Parameter(double val);
  explicit Parameter(const QString &val);

  Parameter(const Parameter &other);
  Parameter &operator=(const Parameter &other);

  Type type() const;
  bool isBool() const;
  bool asBool() const;
  bool isInt() const;
  int asInt() const;
  bool isFloat() const;
  double asFloat() const;
  bool isString() const;
  QString asString() const;

  QDomElement serialize(QDomDocument &doc) const;

public:
  static Parameter fromXml(const QDomElement &node, ParserInfo &info);

protected:
  Type _type;
  QVariant _value;
};


class ParameterEdit: public QWidget
{
  Q_OBJECT

public:
  ParameterEdit(const Parameter &param, QWidget *parent=0);

  Parameter parameter() const;

protected:
  Parameter _parameter;
  QWidget *_edit;
};


class NodeConfigDialog: public QDialog
{
  Q_OBJECT

public:
  NodeConfigDialog(NodeBase *node);

protected slots:
  void apply();

protected:
  NodeBase *_node;
  QLineEdit *_label;
  QTextEdit *_description;
  QHash<QString, ParameterEdit *> _params;
};


#endif // PARAMETER_HH
