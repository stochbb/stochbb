#ifndef NODES_HH
#define NODES_HH

#include "qnetview.hh"
#include <QHash>
#include <QDomElement>
#include <QDialog>
#include "parameter.hh"
#include "lib/api.hh"


class Network;
class Assembler;
class QLineEdit;


class Socket: public QNetSocket
{
  Q_OBJECT

public:
  Socket(QNetSocket::Side side, const QString &name, const QString &label, QNetNode *node);

  const QString &name() const;

protected:
  QString _name;
};


class NodeBase: public QNetNode
{
  Q_OBJECT

public:
  explicit NodeBase(const QString &label, QNetView *parent=0);

  const QString &id() const;
  const QString &type() const;

  bool hasParameters() const;
  bool hasParameter(const QString &name) const;
  Parameter parameter(const QString &name) const;
  const QHash<QString, Parameter> &parameters() const;
  void addParameter(const QString &name, const Parameter &param);
  virtual bool setParameter(const QString &name, const Parameter &param);

  bool hasSocket(const QString &name) const;
  Socket *socket(const QString &name) const;
  void addSocket(QNetSocket *socket);
  virtual QDomElement serialize(QDomDocument &doc) const;

  virtual bool needsPreprocessing(Assembler &assembler) const;
  virtual bool preprocess(Assembler &assembler) const;
  virtual bool processable(Assembler &assembler) const;
  virtual bool assemble(Assembler &assembler) const = 0;

public:
  static NodeBase *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);

protected:
  QString _id;
  QString _type;
  QHash<QString, Socket *> _sockets;
  QHash<QString, Parameter> _params;

  typedef NodeBase *(*nodeFactoryFunction)(const QDomElement &, ParserInfo &, QHash<QString, NodeBase *> &);
  static QHash<QString, nodeFactoryFunction> _factoryFunctions;
};


class TriggerNode: public NodeBase
{
  Q_OBJECT

public:
  TriggerNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static TriggerNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class DelayNode: public NodeBase
{
  Q_OBJECT

public:
  DelayNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static DelayNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class RandomDelayNode: public NodeBase
{
  Q_OBJECT

public:
  RandomDelayNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static RandomDelayNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class GammaProcessNode: public NodeBase
{
  Q_OBJECT

public:
  GammaProcessNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static GammaProcessNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class CompoundGammaProcessNode: public NodeBase
{
  Q_OBJECT

public:
  CompoundGammaProcessNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static CompoundGammaProcessNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class InvGammaProcessNode: public NodeBase
{
  Q_OBJECT

public:
  InvGammaProcessNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static InvGammaProcessNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class CompoundInvGammaProcessNode: public NodeBase
{
  Q_OBJECT

public:
  CompoundInvGammaProcessNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static CompoundInvGammaProcessNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class WeibullProcessNode: public NodeBase
{
  Q_OBJECT

public:
  WeibullProcessNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static WeibullProcessNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class CompoundWeibullProcessNode: public NodeBase
{
  Q_OBJECT

public:
  CompoundWeibullProcessNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static CompoundWeibullProcessNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class MinimumNode: public NodeBase
{
  Q_OBJECT

public:
  MinimumNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static MinimumNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class MaximumNode: public NodeBase
{
  Q_OBJECT

public:
  MaximumNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static MaximumNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class InhibitionNode: public NodeBase
{
  Q_OBJECT

public:
  InhibitionNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  bool assemble(Assembler &assembler) const;

public:
  static InhibitionNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class JoinNode: public NodeBase
{
  Q_OBJECT

public:
  JoinNode(InhibitionNode *sibling, Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;

  bool assemble(Assembler &assembler) const;

public:
  static JoinNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);

protected:
  InhibitionNode *_sibling;
};


class AffineNode: public NodeBase
{
  Q_OBJECT

public:
  AffineNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static AffineNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class ConstantNode: public NodeBase
{
  Q_OBJECT

public:
  ConstantNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static ConstantNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class GammaVarNode: public NodeBase
{
  Q_OBJECT

public:
  GammaVarNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static GammaVarNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class CompoundGammaVarNode: public NodeBase
{
  Q_OBJECT

public:
  CompoundGammaVarNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static CompoundGammaVarNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class InvGammaVarNode: public NodeBase
{
  Q_OBJECT

public:
  InvGammaVarNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static InvGammaVarNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class CompoundInvGammaVarNode: public NodeBase
{
  Q_OBJECT

public:
  CompoundInvGammaVarNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static CompoundInvGammaVarNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class WeibullVarNode: public NodeBase
{
  Q_OBJECT

public:
  WeibullVarNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static WeibullVarNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class CompoundWeibullVarNode: public NodeBase
{
  Q_OBJECT

public:
  CompoundWeibullVarNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static CompoundWeibullVarNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class UniformVarNode: public NodeBase
{
  Q_OBJECT

public:
  UniformVarNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static UniformVarNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class NormalVarNode: public NodeBase
{
  Q_OBJECT

public:
  NormalVarNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static NormalVarNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class CompoundNormalVarNode: public NodeBase
{
  Q_OBJECT

public:
  CompoundNormalVarNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;
  virtual bool assemble(Assembler &assembler) const;

public:
  static CompoundNormalVarNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class OutputNode: public NodeBase
{
  Q_OBJECT

protected:
  OutputNode(const QString &label, QNetView *parent=0);
  virtual bool assemble(Assembler &assembler) const;

public:
  virtual void execute(const QHash<Socket *, stochbb::Var> &vartable) = 0;
};


class MarginalPlotNode: public OutputNode
{
  Q_OBJECT

public:
  MarginalPlotNode(Network *parent=0);

  bool setParameter(const QString &name, const Parameter &param);
  QDomElement serialize(QDomDocument &doc) const;

  void execute(const QHash<Socket *, stochbb::Var> &vartable);

public:
  static MarginalPlotNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class ScatterPlotNode: public OutputNode
{
  Q_OBJECT

public:
  ScatterPlotNode(Network *parent=0);

  QDomElement serialize(QDomDocument &doc) const;

  void execute(const QHash<Socket *, stochbb::Var> &vartable);

public:
  static ScatterPlotNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


class KDEPlotNode: public OutputNode
{
  Q_OBJECT

public:
  KDEPlotNode(Network *parent=0);

  bool setParameter(const QString &name, const Parameter &param);
  QDomElement serialize(QDomDocument &doc) const;

  void execute(const QHash<Socket *, stochbb::Var> &vartable);

public:
  static KDEPlotNode *fromXml(const QDomElement &node, ParserInfo &info, QHash<QString, NodeBase *> &nodeTable);
};


#endif // NODES_HH
