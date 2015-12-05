#include "xmlparser.hh"
#include "chain.hh"
#include "minmax.hh"
#include "mixture.hh"
#include "exception.hh"
#include <iostream>
#include <QFile>
#include <QFileInfo>


using namespace sbb;

/* ********************************************************************************************* *
 * Implementation of VariableDefinition
 * ********************************************************************************************* */
XmlParser::VariableDefinition::VariableDefinition()
{
  // pass...
}

XmlParser::VariableDefinition::~VariableDefinition() {
  // pass...
}


/* ********************************************************************************************* *
 * Implementation of GenericVariableDefinition
 * ********************************************************************************************* */
XmlParser::GenericVariableDefinition::GenericVariableDefinition(
    Var (*func)(QDomElement &, ContextObj *, XmlParser *))
  : VariableDefinition(), _func(func)
{
  // pass...
}

Var
XmlParser::GenericVariableDefinition::instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser) {
  return _func(node, ctx, parser);
}


/* ********************************************************************************************* *
 * Implementation of UserVariableDefinition
 * ********************************************************************************************* */
XmlParser::UserVariableDefinition::UserVariableDefinition(QDomElement &node)
  : VariableDefinition(), _definition(node.cloneNode().toElement()), _name()
{
  // Get name
  if (! _definition.hasAttribute("name")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": cannot create variable defintion: No name specified.";
    throw err;
  }
  _name = _definition.attribute("name");
}

Var
XmlParser::UserVariableDefinition::instantiate(QDomElement &node, ContextObj *ctx, XmlParser *parser) {
  // Get parameter specification
  QHash<QString, double> params = parser->parseParams(node, ctx);
  // Create local context of variable definition
  ContextObj *local = new ContextObj(ctx);
  /* Process definition */
  for (QDomElement p=_definition.firstChildElement(); !p.isNull(); p=p.nextSiblingElement()) {
    if ("param" == p.tagName()) {
      if (! p.hasAttribute("id")) {
        ParserError err;
        err << "ParserError @" << p.lineNumber()
            << ": <param> has no 'id' attribute.";
        throw err;
      }
      QString id = p.attribute("id");
      if (! params.contains(id)) {
        ParserError err;
        err << "ParserError @" << p.lineNumber()
            << ": Cannot instantiate variable '" << _name.toStdString()
            << "', parameter '" << id.toStdString() << "' not given.";
        throw err;
      }
      // Store parameter locally such that it is available to the expressions
      // for the variable instantiations
      local->setParam(id.toStdString(), params[id]);
    } else if ("var" == p.tagName()) {
      parser->parseVarDef(p, local);
    } else {
      ParserError err;
      err << "ParserError @" << p.lineNumber()
          << ": unexpected tag '" << p.tagName().toStdString() << "'. Expected 'param' or 'var'.";
      throw err;
    }
  }

  // Check if there exists a variable in local context which matches the definition name
  if (! local->hasVar(_name.toStdString())) {
    ParserError err;
    err << "ParserError: Cannot instantiate variable " << _name.toStdString()
        << ": No variable named " << _name.toStdString() << " defined.";
    throw err;
  }

  Var var(local->var(_name.toStdString()));
  if (node.hasAttribute("name")) { var.setName(node.attribute("name").toStdString()); }
  else if (node.hasAttribute("id")) { var.setName(node.attribute("id").toStdString()); }
  return var;
}


/* ********************************************************************************************* *
 * Implementation of XmlParser
 * ********************************************************************************************* */
XmlParser::XmlParser() {
  _factories["delta"] = new GenericVariableDefinition(&XmlParser::parseDelta);
  _factories["unif"] = new GenericVariableDefinition(&XmlParser::parseUnif);
  _factories["normal"] = new GenericVariableDefinition(&XmlParser::parseNorm);
  _factories["gamma"] = new GenericVariableDefinition(&XmlParser::parseGamma);
  _factories["chain"] = new GenericVariableDefinition(&XmlParser::parseChain);
  _factories["maximum"] = new GenericVariableDefinition(&XmlParser::parseMaximum);
  _factories["minimum"] = new GenericVariableDefinition(&XmlParser::parseMinimum);
  _factories["mixture"] = new GenericVariableDefinition(&XmlParser::parseMixture);
}

XmlParser::~XmlParser() {
  QHash<QString, VariableDefinition *>::iterator item = _factories.begin();
  for (; item != _factories.end(); item++) {
    delete item.value();
  }
}

Simulation
XmlParser::parse(const QString &filename) {
  QFile file(filename);
  if (! file.open(QIODevice::ReadOnly)) {
    ParserError err;
    err << "Cannot open file " << filename.toStdString() << std::endl;
    throw err;
  }

  QDomDocument doc;
  QString msg; int row;
  if (! doc.setContent(&file, true, &msg, &row)) {
    file.close();
    ParserError err;
    err << "Cannot parse file " << filename.toStdString() << ": " << msg.toStdString()
        << " @" << row << std::endl;
    throw err;
  }

  // Store file-path as current path
  QFileInfo info(file.fileName());
  _pathStack.push_back(info.absolutePath());
  QDomElement root = doc.documentElement();
  return parse(root);
}

Simulation
XmlParser::parse(QDomElement &root) {
  SimulationObj *sim = new SimulationObj();
  if ("simulation" != root.tagName()) { return sim; }
  // Parse elements
  for (QDomElement node=root.firstChildElement(); !node.isNull(); node=node.nextSiblingElement()) {
    // dispatch by element name
    if ("var" == node.tagName()) { parseVarDef(node, sim); }
    else if ("output" == node.tagName()) { parseOutput(node, sim); }
    else if ("define" == node.tagName()) { parseDefine(node, sim); }
    else if ("param" == node.tagName())  { parseSimParam(node, sim); }
    else if ("load" == node.tagName()) { parseLoad(node, sim); }
    else {
      ParserError err;
      err << "ParserError @"<< node.lineNumber()
          << ": Unexpected element " << node.tagName().toStdString();
      throw err;
    }
  }
  return sim;
}

void
XmlParser::parseLoad(QDomElement &node, ContextObj *ctx) {
  QString name = node.text();
  QFile file(_pathStack.back() + "/" + name + ".xml");
  if (! file.open(QIODevice::ReadOnly)) {
    ParserError err;
    err << "Cannot import file " << file.fileName().toStdString() << ": File not readable.";
    throw err;
  }

  QString message; int line;
  QDomDocument doc;
  if (! doc.setContent(&file, true, &message, &line)) {
    ParserError err;
    err << "ParserError: " << file.fileName().toStdString() << "@" << line
        << ": " << message.toStdString();;
    throw err;
  }

  QDomElement root = doc.documentElement();
  if ("module" != root.tagName()) {
    ParserError err;
    err << "ParserError: " << file.fileName().toStdString() << "@" << root.lineNumber()
        << ": Unexpected tag <" << root.tagName().toStdString() << ">. Expected <module>";
    throw err;
  }

  QFileInfo info(file.fileName());
  _pathStack.push_back(info.absolutePath());
  try { parseModule(root, ctx); }
  catch (...) { _pathStack.pop_back(); throw;}
  _pathStack.pop_back();
}

void
XmlParser::parseModule(QDomElement &node, ContextObj *ctx) {
  for (QDomElement p=node.firstChildElement(); !p.isNull(); p=p.nextSiblingElement()) {
    // Dispatch by type
    if ("define" == p.tagName()) { parseDefine(p, ctx); }
    else if ("load" == p.tagName()) { parseLoad(p, ctx); }
    else {
      ParserError err;
      err << "ParseError @" << p.lineNumber()
          << ": Unexpected tag <" << p.tagName().toStdString()
          << ">, expected <load> or <define>.";
      throw err;
    }
  }
}


Var
XmlParser::parseDelta(QDomElement &node, ContextObj *sim, XmlParser *parser) {
  QHash<QString, double> params = parser->parseParams(node, sim);
  if (0 == params.count("delay")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'delay' parameter.";
    throw err;
  }

  std::string name = "";
  if (node.hasAttribute("name")) { name = node.attribute("name").toStdString(); }
  else if (node.hasAttribute("id")) { name = node.attribute("id").toStdString(); }
  return AtomicVarObj::delta(params["delay"], name);
}

Var
XmlParser::parseUnif(QDomElement &node, ContextObj *sim, XmlParser *parser) {
  QHash<QString, double> params = parser->parseParams(node, sim);
  if (0 == params.count("a")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'a' parameter.";
    throw err;
  }
  if (0 == params.count("b")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'b' parameter.";
    throw err;
  }

  std::string name = "";
  if (node.hasAttribute("name")) { name = node.attribute("name").toStdString(); }
  else if (node.hasAttribute("id")) { name = node.attribute("id").toStdString(); }
  return AtomicVarObj::unif(params["a"], params["b"], name);
}

Var
XmlParser::parseNorm(QDomElement &node, ContextObj *sim, XmlParser *parser) {
  QHash<QString, double> params = parser->parseParams(node, sim);
  if (0 == params.count("mu")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'mu' parameter.";
    throw err;
  }
  if (0 == params.count("sigma")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'sigma' parameter.";
    throw err;
  }

  std::string name = "";
  if (node.hasAttribute("name")) { name = node.attribute("name").toStdString(); }
  else if (node.hasAttribute("id")) { name = node.attribute("id").toStdString(); }
  return AtomicVarObj::norm(params["mu"], params["sigma"], name);
}

Var
XmlParser::parseGamma(QDomElement &node, ContextObj *sim, XmlParser *parser) {
  QHash<QString, double> params = parser->parseParams(node, sim);
  if (0 == params.count("k")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'k' parameter.";
    throw err;
  }
  if (0 == params.count("theta")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'theta' parameter.";
    throw err;
  }

  std::string name = "";
  if (node.hasAttribute("name")) { name = node.attribute("name").toStdString(); }
  else if (node.hasAttribute("id")) { name = node.attribute("id").toStdString(); }
  return AtomicVarObj::gamma(params["k"], params["theta"], name);
}

Var
XmlParser::parseChain(QDomElement &node, ContextObj *sim, XmlParser *parser) {
  QVector<Var> vars = parser->parseVars(node, sim);
  if (0 == vars.size()) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no variables.";
    throw err;
  }

  std::string name = "";
  if (node.hasAttribute("name")) { name = node.attribute("name").toStdString(); }
  else if (node.hasAttribute("id")) { name = node.attribute("id").toStdString(); }
  return new ChainObj(vars.toStdVector(), name);
}

Var
XmlParser::parseMaximum(QDomElement &node, ContextObj *sim, XmlParser *parser) {
  QVector<Var> vars = parser->parseVars(node, sim);
  if (0 == vars.size()) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no variables.";
    throw err;
  }

  std::string name = "";
  if (node.hasAttribute("name")) { name = node.attribute("name").toStdString(); }
  else if (node.hasAttribute("id")) { name = node.attribute("id").toStdString(); }
  return new MaximumObj(vars.toStdVector(), name);
}

Var
XmlParser::parseMinimum(QDomElement &node, ContextObj *sim, XmlParser *parser) {
  QVector<Var> vars = parser->parseVars(node, sim);
  if (0 == vars.size()) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no variables.";
    throw err;
  }

  std::string name = "";
  if (node.hasAttribute("name")) { name = node.attribute("name").toStdString(); }
  else if (node.hasAttribute("id")) { name = node.attribute("id").toStdString(); }
  return new MinimumObj(vars.toStdVector(), name);
}

Var
XmlParser::parseMixture(QDomElement &node, ContextObj *sim, XmlParser *parser) {
  std::vector< std::pair<double, Var> > vars = parser->parseWeights(node, sim);
  if (0 == vars.size()) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no variables.";
    throw err;
  }

  std::string name = "";
  if (node.hasAttribute("name")) { name = node.attribute("name").toStdString(); }
  else if (node.hasAttribute("id")) { name = node.attribute("id").toStdString(); }
  // Get weights and variables
  std::vector<double> weights; weights.reserve(vars.size());
  std::vector<Var> variables; variables.reserve(vars.size());
  for (size_t i=0; i<vars.size(); i++) {
    weights.push_back(vars[i].first);
    variables.push_back(vars[i].second);
  }
  // done
  return new MixtureObj(weights, variables, name);
}

Var
XmlParser::parseVar(QDomElement &node, ContextObj *sim) {
  if ("var" == node.tagName()) {
    if (node.hasAttribute("type") && (! node.hasAttribute("ref"))) {
      return parseVarDef(node, sim);
    } else if ( (! node.hasAttribute("type")) && node.hasAttribute("ref")) {
      return parseVarRef(node, sim);
    }
    ParserError err;
    err << "ParserError @"<< node.lineNumber()
        << ": <var> element must have either a type or a ref attribute.";
    throw err;
  }

  ParserError err;
  err << "ParserError @"<< node.lineNumber()
      << ": Unexpected element " << node.tagName().toStdString()
      << " expected <var>.";
  throw err;
}

Var
XmlParser::parseVarDef(QDomElement &node, ContextObj *sim) {
  if (! node.hasAttribute("type")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no type attribute.";
    throw err;
  }
  if (! _factories.contains(node.attribute("type"))) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << ": Unknown variable type '"
        << node.attribute("type").toStdString() << "'.";
    throw err;
  }
  Var var = _factories[node.attribute("type")]->instantiate(node, sim, this);
  if (node.hasAttribute("id")) {
    sim->addVar(node.attribute("id").toStdString(), var);
  }
  return var;
}

Var
XmlParser::parseVarRef(QDomElement &node, ContextObj *sim) {
  QString name = node.attribute("ref");
  if (! sim->hasVar(name.toStdString())) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString()
        << ": variable " << name.toStdString() << " not known.";
    throw err;
  }
  return sim->var(name.toStdString());
}

QHash<QString, double>
XmlParser::parseParams(QDomElement &node, ContextObj *sim) {
  QHash<QString, double> params;
  for (QDomElement p=node.firstChildElement("param"); !p.isNull(); p=p.nextSiblingElement("param")) {
    if (! p.hasAttribute("name")) {
      ParserError err;
      err << "ParserError @" << p.lineNumber()
          << ": " << node.tagName().toStdString() << " has no 'name' attribute.";
      throw err;
    }
    QDomElement mml = p.firstChildElement();
    params.insert(p.attribute("name"), parseMathML(mml, sim));
  }
  return params;
}

QVector<Var> XmlParser::parseVars(QDomElement &node, ContextObj *sim) {
  QVector<Var> vars;
  for (QDomElement p=node.firstChildElement(); !p.isNull(); p=p.nextSiblingElement()) {
    vars.push_back(parseVar(p, sim));
  }
  return vars;
}

std::vector< std::pair<double, Var> >
XmlParser::parseWeights(QDomElement &node, ContextObj *sim) {
  std::vector< std::pair<double, Var> > weights;
  for (QDomElement p=node.firstChildElement(); !p.isNull(); p=p.nextSiblingElement()) {
    weights.push_back(parseWeight(p, sim));
  }
  return weights;
}

std::pair<double, Var>
XmlParser::parseWeight(QDomElement &node, ContextObj *sim) {
  if ("weight" != node.tagName()) {
    ParserError err;
    err << "ParserError @"<< node.lineNumber()
        << ": Unexpected element <" << node.tagName().toStdString()
        << ">. Expected <weight>.";
    throw err;
  }
  // get weight
  QDomElement wnode = node.firstChildElement();
  double w = parseMathML(wnode, sim);
  // get variable
  QDomElement vnode = wnode.nextSiblingElement();
  Var var = parseVar(vnode, sim);
  return std::pair<double, Var>(w, var);
}

void
XmlParser::parseDefine(QDomElement &node, ContextObj *ctx) {
  if (! node.hasAttribute("name")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'name' attribute.";
    throw err;
  }
  QString name = node.attribute("name");
  if(_factories.contains(name)) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": cannot redefine variable type '"<< name.toStdString() <<"'.";
    throw err;
  }
  _factories[name] = new UserVariableDefinition(node);
}

void
XmlParser::parseOutput(QDomElement &node, SimulationObj *sim) {
  if (! node.hasAttribute("from")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'from' attribute.";
    throw err;
  }
  if (! node.hasAttribute("to")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'to' attribute.";
    throw err;
  }
  if (! node.hasAttribute("steps")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'steps' attribute.";
    throw err;
  }
  sim->setTMin(node.attribute("from").toDouble());
  sim->setTMax(node.attribute("to").toDouble());
  sim->setSteps(node.attribute("steps").toUInt());
  // Parse output variables
  for (QDomElement p=node.firstChildElement(); !p.isNull(); p=p.nextSiblingElement()) {
    sim->addOutputVar(parseVar(p, sim));
  }
}

void
XmlParser::parseSimParam(const QDomElement &node, SimulationObj *sim) {
  if (! node.hasAttribute("id")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'id' attribute.";
    throw err;
  }
  sim->setParam(node.attribute("id").toStdString(),
                parseMathML(node.firstChildElement(), sim));
}

double
XmlParser::parseMathML(const QDomElement &node, ContextObj *ctx) {
  // Dispatch by tag name
  if ("cn" == node.tagName()) {
    return parseMMLNumber(node, ctx);
  } else if ("ci" == node.tagName()) {
    return parseMMLSymbol(node, ctx);
  } else if ("apply" == node.tagName()) {
    return parseMMLApply(node, ctx);
  }
  ParserError err;
  err << "ParserError @" << node.lineNumber()
      << " : Unknown MathML element '" << node.tagName().toStdString()
      << ". Expected <ci>, <cn> or <apply>.";
  throw err;
}

double
XmlParser::parseMMLNumber(const QDomElement &node, ContextObj *ctx) {
  bool ok=true;
  double value = node.text().toDouble(&ok);
  if (!ok) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": Value '" << node.text().toStdString() << "' is not a number!";
    throw err;
  }
  return value;
}

double
XmlParser::parseMMLSymbol(const QDomElement &node, ContextObj *ctx) {
  std::string name = node.text().toStdString();
  if (! ctx->hasParam(name)) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": Symbol (parameter) '" << name << "' not defined.";
    throw err;
  }
  return ctx->param(name);
}

double
XmlParser::parseMMLApply(const QDomElement &node, ContextObj *ctx) {
  if (! node.hasChildNodes()) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": Empty <apply> tag.";
    throw err;
  }
  // Get operator
  QDomElement op = node.firstChildElement();
  std::vector<double> args; args.reserve(10);
  // Collect arguments
  for (QDomElement p = op.nextSiblingElement(); ! p.isNull(); p=p.nextSiblingElement()) {
    args.push_back(parseMathML(p, ctx));
  }
  // dispatch by operator type
  if ("plus" == op.tagName()) {
    double res = 0;
    for (size_t i=0; i<args.size(); i++) {
      res += args[i];
    }
    return res;
  } else if ("times" == op.tagName()) {
    double res = 1;
    for (size_t i=0; i<args.size(); i++) {
      res *= args[i];
    }
    return res;
  } else if ("minus" == op.tagName()) {
    if (1 == args.size()) {
      return -args[0];
    } else if (2 != args.size()) {
      return args[0]-args[1];
    }
    ParserError err;
    err << "ParserError @" << op.lineNumber()
        << ": the <minus> operator requires exactly 1 or 2 arguments.";
    throw err;
  } else if ("divide" == op.tagName()) {
    if (2 != args.size()) {
      ParserError err;
      err << "ParserError @" << op.lineNumber()
          << ": the <divide> operator requires exactly 2 arguments.";
      throw err;
    }
    return args[0]/args[1];
  } else if ("power" == op.tagName()) {
    if (2 != args.size()) {
      ParserError err;
      err << "ParserError @" << op.lineNumber()
          << ": the <power> operator requires exactly 2 arguments.";
      throw err;
    }
    return std::pow(args[0],args[1]);
  } else if ("log" == op.tagName()) {
    if (2 != args.size()) {
      ParserError err;
      err << "ParserError @" << op.lineNumber()
          << ": the <log> function requires exactly 1 argument.";
      throw err;
    }
    return std::log(args[0]);
  } else if ("exp" == op.tagName()) {
    if (2 != args.size()) {
      ParserError err;
      err << "ParserError @" << op.lineNumber()
          << ": the <exp> function requires exactly 1 argument.";
      throw err;
    }
    return std::exp(args[0]);
  }

  ParserError err;
  err << "ParserError @" << op.lineNumber()
      << ": Unknown operator or function <" << op.tagName().toStdString() << ">.";
  throw err;
}
