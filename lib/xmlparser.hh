/** \page xml Process & simulation description in XML
 * The @c XmlParser class implements a paser that allows to specifiy process descriptions in XML.
 *
 * A simulation is a simple collection of random variables together with the specification which
 * variables as used for the output. For example the XML code
 * \code
 *  <?xml version="1.0"?>
 *  <simulation xmlns="http://hmatuschek.github.io/stochbb/simulation-0.0.dtd"
 *              xmlns:m="http://www.w3.org/1998/Math/MathML">
 *
 *   <var type="gamma" id="gamma">
 *    <param name="k"> <m:cn>5</m:cn> </param>
 *    <param name="theta"> <m:cn>30</m:cn> </param>
 *   </var>
 *
 *   <output from="0" to="100">
 *    <var ref="gamma"/>
 *   </output>
 *
 *  </simulation>
 * \endcode
 * specifies a trivial simulation.
 *
 * Here a single gamma-distributed random variable is defined.
 * It gets the identifier "gamma". With this identifier, it is possible to refer to the random
 * variable later. The @c type attribute specifies the type of the random variable. Here the
 * build-in type "gamma" is used, specifying a gamma-distributed random variable. The shape and
 * scale parameters, "k" and "theta", are specified using the "param" elements. The "name" attribute
 * of each param element specified the name of the parameter while its value is given as an MathML
 * child element of the param element.
 *
 * Finally, the "output" element specifies which random variables are used for output. Each child
 * element of the output element specifies one variable. Either by referencing a previously defined
 * random variable or by specifying the random variable directly. For example,
 * \code
 *  <?xml version="1.0"?>
 *  <simulation xmlns="http://hmatuschek.github.io/stochbb/simulation-0.0.dtd"
 *              xmlns:m="http://www.w3.org/1998/Math/MathML">
 *
 *   <output from="0" to="100">
 *    <var type="gamma">
 *     <param name="k"> <m:cn>5</m:cn> </param>
 *     <param name="theta"> <m:cn>30</m:cn> </param>
 *    </var>
 *   </output>
 *
 *  </simulation>
 * \endcode
 * is fully equivalent to the previous example.
 *
 * \section xmlrvs Build-in random variables
 * There are only very few build-in random variable types. They can be devided into two groups:
 * atomic and derived random variables.
 *
 * \subsection xmlatomic Atomic random variables
 * Atomic random variables are RVs which follow a specific
 * distribution and are independent from any other RV defined.
 *
 *   Type | Parameters | Process description
 *   --- | --- | ---
 *   @c delta | delay | A constant delay or a process with a fixed waiting time.
 *   @c uniform | a, b | A process with a uniform-distributed waiting time.
 *   @c normal | mu, sigma | A process with a normal-distributed waiting.
 *   @c gamma | k, theta | A process with a gamma-distributed waiting time.
 *
 * Beside these atomic random variables, there are derived random variables which depend on other
 * RVs. Hence they usually to not taks simple parameters as arguments but other random variables.
 *
 * \subsection xmlchain Chains of random variables.
 * One of the most important derived random variable is a chain. A chain can be defined as
 * \code
 *  <var type="chain">
 *   <var ref="X1"/> <var ref="X2"/> <var ref="X2"/>
 *  </var>
 * \endcode
 * Here the random variable X1, X2 and X3 are chained. This means, for example, that X1 triggers
 * X2 and X2 then triggers the process X3. Mathematically, the waiting time of the chain is simply
 * the sum of the waiting times of the processes X1-X3.
 *
 * \subsection xmlmax Maximum of parallel processes
 * Another derived random variable is the "maximum" type. This type can be interpreted as a
 * processes that waits until all referred parallel processes are completed hence the waiting time
 * of
 * \code
 *  <var type="maximum">
 *   <var ref="X1"/> <var ref="X2"/> <var ref="X2"/>
 *  </var>
 * \endcode
 * is just the maximum of the random variables X1-X3.
 *
 * \subsection xmlmin Minimum of parallel processes
 * Likewise, the "minimum" derived random variable is simply the minimum of the referred random
 * variables and represents the waiting time of the fastest parallel processes X1-X3:
 * \code
 *  <var type="minimum">
 *   <var ref="X1"/> <var ref="X2"/> <var ref="X2"/>
 *  </var>
 * \endcode
 *
 * \section xmluserrv User defined random variable types
 * The limited build-in random variable type would turn the definition of complex stationary random
 * processess difficult. Hence it is possible to define new typed derived from the build-in ones.
 * For example
 * \code
 *  <define id="exp">
 *   <param id="lambda"/>
 *
 *   <var type="gamma" id="exp">
 *    <param name="k"><m:cn>1</m:cn></param>
 *    <param name="theta">
 *     <m:apply>
 *      <m:divide/>
 *      <m:cn>1</m:cn>
 *      <m:ci>lambda</m:ci>
 *     </m:apply>
 *    </param>
 *   </var>
 *  </define>
 * \endcode
 * defines a new type "exp" as a special case of the Gamma distribution as
 * \f$Exp(\lambda) = \Gamma(1,\lambda^{-1}\f$. The define takes an "id" attribute specifying the
 * new type ID. Then one ore more "param" elments are used to define the paramters of the new type.
 * Here only one parameter named "lambda" is defined. Finally a random variable with the same
 * identifier as the type name is defined which implements the derived random variable type.
 *
 * A more complex example might by
 * \code
 *  <define id="exgauss">
 *   <param id="lambda">
 *   <param id="mu"/>
 *   <param id="sigma"/>
 *
 *   <var type="gamma" id="ex">
 *    <param name="k"><m:cn>1</m:cn></param>
 *    <param name="theta">
 *     <m:apply>
 *      <m:divide/>
 *      <m:cn>1</m:cn>
 *      <m:ci>lambda</m:ci>
 *     </m:apply>
 *    </param>
 *   </var>
 *
 *   <var type="normal" id="gauss">
 *    <param name="mu"> <m:ci>mu</m:ci> </param>
 *    <param name="sigma"> <m:ci>sigma</m:ci> </param>
 *   </var>
 *
 *   <var type="chain" id="exgauss">
 *     <var ref="ex"/> <var ref="gauss"/>
 *   </var>
 *  </define>
 * \endcode
 * This defines a new random variable type with an Ex-Gauss distribution by chaining a
 * exponential (defined as a specialcase of the Gamma distribution) and a guassian distributed
 * random process. This new type can then be used to instantiate Ex-Gauss distributed random
 * varaibles.
 * \code
 *  <var type="exgauss">
 *   <param name="lambda"> <m:cn>10</m:cn> </param>
 *   <param name="mu"> <m:cn>20</m:cn> </param>
 *   <param name="sigma"> <m:cn>10</m:cn> </param>
 *  </var>
 * \endcode
 *
 * \section xmlmod Modules
 * Some of the aforementioned user defined random variable types may be reused in several
 * simulations. For these cases it is convenient to export these definitions into separate files
 * and load them when needed. These files are called modules and form simple collections of
 * variable type definitions. For example a module could look like
 * \code
 *  <?xml version="1.0"?>
 *  <module xmlns="http://hmatuschek.github.io/stochbb/module-0.0.dtd"
 *              xmlns:m="http://www.w3.org/1998/Math/MathML">
 *
 *  <define id="type">
 *   ...
 *  </define>
 *  </module>
 * \endcode
 * Within the body of the module, the user defined types can be spcified. Such a module can then
 * be loaded using the "load" tag within a simulation. Please note that it is also possible to
 * load other modules from within a module.
 * \code
 *  <simulation xmlns="http://hmatuschek.github.io/stochbb/simulation-0.0.dtd"
 *              xmlns:m="http://www.w3.org/1998/Math/MathML">
 *
 *   <load>FILENAME_W/O_EXTENSION</load>
 *
 *   <output ...>
 *     ...
 *   </output>
 *  </simulation>
 * \endcode
 */

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
  /** Parses a 'param' statement whithin a simulation. */
  void parseSimParam(const QDomElement &node, SimulationObj *sim);

  /** Parses a list of variable reference statements. */
  QVector<Var> parseVars(QDomElement &node, ContextObj *sim);
  /** Parses a variable reference, either a 'ref' or 'var' statement. */
  Var parseVar(QDomElement &node, ContextObj *sim);
  /** Parses a 'var' statement. */
  Var parseVarDef(QDomElement &node, ContextObj *sim);
  /** Parses a 'ref' statement. */
  Var parseVarRef(QDomElement &node, ContextObj *sim);
  /** Parses a list of 'param' statements. */
  QHash<QString, double> parseParams(QDomElement &node, ContextObj *sim);

  /** Parses and evaluates MathML expressions. */
  double parseMathML(const QDomElement &node, ContextObj *ctx);
  /** Parses a MathML 'cn' expression. */
  double parseMMLNumber(const QDomElement &node, ContextObj *ctx);
  /** Parses a MathML 'ci' expression. */
  double parseMMLSymbol(const QDomElement &node, ContextObj *ctx);
  /** Parses a MathML 'apply' expression. */
  double parseMMLApply(const QDomElement &node, ContextObj *ctx);

  /** Instantiates a delta distribution. */
  static Var parseDelta(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a uniform distribution. */
  static Var parseUnif(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a normal distribution. */
  static Var parseNorm(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a gamma distribution. */
  static Var parseGamma(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a chain of random processes. */
  static Var parseChain(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a maximum of random variables. */
  static Var parseMaximum(QDomElement &node, ContextObj *sim, XmlParser *parser);
  /** Instantiates a minimum of random variables. */
  static Var parseMinimum(QDomElement &node, ContextObj *sim, XmlParser *parser);

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
    virtual Var instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser) = 0;
  };

  /** User defined random variable definition in XML.
   * @ingroup internal */
  class UserVariableDefinition: public VariableDefinition
  {
  public:
    /** Constructor. */
    UserVariableDefinition(QDomElement &node);

    virtual Var instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser);

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
    GenericVariableDefinition(Var (*func)(QDomElement &node, ContextObj *ctx, XmlParser *parser));

    Var instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser);

  protected:
    /** The factory function. */
    Var (*_func)(QDomElement &node, ContextObj *ctx, XmlParser *parser);
  };

protected:
  /** Table of known variable types. */
  QHash<QString, VariableDefinition *> _factories;
  /** Stack of current directories, needed to process 'load' statements with relative paths. */
  QList<QString> _pathStack;
};

}

#endif // __SBB_XMLPARSER_HH__
