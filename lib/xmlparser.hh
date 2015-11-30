#ifndef __SBB_XMLPARSER_HH__
#define __SBB_XMLPARSER_HH__

#include "api.hh"

#include <QDomElement>
#include <QVector>
#include <QHash>

namespace sbb {

// Forward decl.
class XmlParser;

class VariableDefinition
{
protected:
  VariableDefinition();

public:
  virtual ~VariableDefinition();
  virtual VarObj *instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser) = 0;
};


class UserVariableDefinition: public VariableDefinition
{
public:
  UserVariableDefinition(QDomElement &node);

  virtual VarObj *instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser);

protected:
  QDomElement _definition;
  QString _name;
};


class GenericVariableDefinition: public VariableDefinition
{
public:
  GenericVariableDefinition(VarObj *(*func)(QDomElement &node, ContextObj *ctx, XmlParser *parser));

  VarObj *instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser);

protected:
  VarObj *(*_func)(QDomElement &node, ContextObj *ctx, XmlParser *parser);
};


class XmlParser
{
public:
  XmlParser();
  ~XmlParser();
  Simulation parse(const QString &filename);
  Simulation parse(QDomElement &doc);

  void parseLoad(QDomElement &node, ContextObj *ctx);
  void parseModule(QDomElement &node, ContextObj *ctx);
  void parseDefine(QDomElement &node, ContextObj *ctx);

  VarObj *parseVar(QDomElement &node, ContextObj *sim);
  VarObj *parseVarDef(QDomElement &node, ContextObj *sim);
  VarObj *parseVarRef(QDomElement &node, ContextObj *sim);

  QHash<QString, double> parseParams(QDomElement &node, ContextObj *sim);
  double parseMathML(QDomElement &node, ContextObj *ctx);

protected:
  void parseOutput(QDomElement &node, SimulationObj *sim);
  QVector<VarObj *> parseVars(QDomElement &node, ContextObj *sim);

  double parseMMLNumber(QDomElement &node, ContextObj *ctx);
  double parseMMLSymbol(QDomElement &node, ContextObj *ctx);
  double parseMMLApply(QDomElement &node, ContextObj *ctx);

  static VarObj *parseDelta(QDomElement &node, ContextObj *sim, XmlParser *parser);
  static VarObj *parseUnif(QDomElement &node, ContextObj *sim, XmlParser *parser);
  static VarObj *parseNorm(QDomElement &node, ContextObj *sim, XmlParser *parser);
  static VarObj *parseGamma(QDomElement &node, ContextObj *sim, XmlParser *parser);
  static VarObj *parseChain(QDomElement &node, ContextObj *sim, XmlParser *parser);
  static VarObj *parseMaximum(QDomElement &node, ContextObj *sim, XmlParser *parser);
  static VarObj *parseMinimum(QDomElement &node, ContextObj *sim, XmlParser *parser);

protected:
  QHash<QString, VariableDefinition *> _factories;
  QList<QString> _pathStack;
};

}

#endif // __SBB_XMLPARSER_HH__
