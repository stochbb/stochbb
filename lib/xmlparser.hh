#ifndef __SBB_XMLPARSER_HH__
#define __SBB_XMLPARSER_HH__

#include "api.hh"
#include "simulation.hh"

#include <QDomElement>
#include <QVector>
#include <QHash>

namespace sbb {

// Forward decl.
class XmlParser;

/** Parses simulations specified in XML.
 * @ingroup internal */
class XmlParser
{
public:
  /** Constructor. */
  XmlParser();
  /** Destructor. */
  ~XmlParser();
  /** Parses the given file. */
  Simulation parse(const QString &filename);
  /** Parses the given XML document. */
  Simulation parse(QDomElement &doc);

protected:
  /** Parses a 'load' statement. */
  void parseLoad(QDomElement &node, ContextObj *ctx);
  /** Parses a 'module' statement. */
  void parseModule(QDomElement &node, ContextObj *ctx);
  /** Parses a 'define' statement. */
  void parseDefine(QDomElement &node, ContextObj *ctx);
  /** Parses a 'output' statement. */
  void parseOutput(QDomElement &node, SimulationObj *sim);
  /** Parses a list of variable reference statements. */
  QVector<VarObj *> parseVars(QDomElement &node, ContextObj *sim);
  /** Parses a variable reference, either a 'ref' or 'var' statement. */
  VarObj *parseVar(QDomElement &node, ContextObj *sim);
  /** Parses a 'var' statement. */
  VarObj *parseVarDef(QDomElement &node, ContextObj *sim);
  /** Parses a 'ref' statement. */
  VarObj *parseVarRef(QDomElement &node, ContextObj *sim);
  /** Parses a list of 'param' statements. */
  QHash<QString, double> parseParams(QDomElement &node, ContextObj *sim);
  /** Parses and evaluates MathML expressions. */
  double parseMathML(QDomElement &node, ContextObj *ctx);
  /** Parses a MathML 'cn' expression. */
  double parseMMLNumber(QDomElement &node, ContextObj *ctx);
  /** Parses a MathML 'ci' expression. */
  double parseMMLSymbol(QDomElement &node, ContextObj *ctx);
  /** Parses a MathML 'apply' expression. */
  double parseMMLApply(QDomElement &node, ContextObj *ctx);
  /** Instantiates a delta distribution. */
  static VarObj *parseDelta(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a uniform distribution. */
  static VarObj *parseUnif(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a normal distribution. */
  static VarObj *parseNorm(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a gamma distribution. */
  static VarObj *parseGamma(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a chain of random processes. */
  static VarObj *parseChain(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a maximum of random variables. */
  static VarObj *parseMaximum(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a minimum of random variables. */
  static VarObj *parseMinimum(QDomElement &node, ContextObj *sim, XmlParser *parser);

protected:
  /** Abstract class of known random variable types.
   * @ingroup internal */
  class VariableDefinition
  {
  protected:
    /** Hidden constructor. */
    VariableDefinition();

  public:
    /** Destructor. */
    virtual ~VariableDefinition();

    /** Needs to be implemented to instantiate a random variable. */
    virtual VarObj *instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser) = 0;
  };

  /** User defined random variable definition in XML.
   * @ingroup internal */
  class UserVariableDefinition: public VariableDefinition
  {
  public:
    /** Constructor. */
    UserVariableDefinition(QDomElement &node);

    virtual VarObj *instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser);

  protected:
    /** The variable type definition. */
    QDomElement _definition;
    /** The name of the variable type. */
    QString _name;
  };

  /** Build-in random variable types.
   * @ingroup internal */
  class GenericVariableDefinition: public VariableDefinition
  {
  public:
    /** Constructor. Taks a factor function. */
    GenericVariableDefinition(VarObj *(*func)(QDomElement &node, ContextObj *ctx, XmlParser *parser));

    VarObj *instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser);

  protected:
    /** The factory function. */
    VarObj *(*_func)(QDomElement &node, ContextObj *ctx, XmlParser *parser);
  };

protected:
  /** Table of known variable types. */
  QHash<QString, VariableDefinition *> _factories;
  /** Stack of current directories, needed to process 'load' statements with relative paths. */
  QList<QString> _pathStack;
};

}

#endif // __SBB_XMLPARSER_HH__
