#include "xmlparser.hh"
#include "chain.hh"
#include "minmax.hh"
#include "exception.hh"
#include <iostream>


using namespace sbb;

XmlParser::XmlParser() {
  _factories["delta"] = &XmlParser::parseDelta;
  _factories["unif"] = &XmlParser::parseUnif;
  _factories["normal"] = &XmlParser::parseNorm;
  _factories["gamma"] = &XmlParser::parseGamma;
  _factories["chain"] = &XmlParser::parseChain;
  _factories["maximum"] = &XmlParser::parseMaximum;
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
    else {
      ParserError err;
      err << "ParserError @"<< node.lineNumber()
          << ": Unexpected element " << node.tagName().toStdString();
      throw err;
    }
  }
  return sim;
}

VarObj *
XmlParser::parseDelta(QDomElement &node, SimulationObj *sim, XmlParser *parser) {
  QHash<QString, double> params = parser->parseParams(node, sim);
  if (0 == params.count("delay")) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no 'delay' parameter.";
    throw err;
  }
  return GenericVarObj::delta(params["delay"]);
}

VarObj *
XmlParser::parseUnif(QDomElement &node, SimulationObj *sim, XmlParser *parser) {
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

  return GenericVarObj::unif(params["a"], params["b"]);
}

VarObj *
XmlParser::parseNorm(QDomElement &node, SimulationObj *sim, XmlParser *parser) {
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

  return GenericVarObj::norm(params["mu"], params["sigma"]);
}

VarObj *
XmlParser::parseGamma(QDomElement &node, SimulationObj *sim, XmlParser *parser) {
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
  return GenericVarObj::norm(params["k"], params["theta"]);
}

VarObj *
XmlParser::parseChain(QDomElement &node, SimulationObj *sim, XmlParser *parser) {
  QVector<VarObj *> vars = parser->parseVars(node, sim);
  if (0 == vars.size()) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no variables.";
    throw err;
  }
  return new ChainObj(vars.toStdVector());
}

VarObj *
XmlParser::parseMaximum(QDomElement &node, SimulationObj *sim, XmlParser *parser) {
  QVector<VarObj *> vars = parser->parseVars(node, sim);
  if (0 == vars.size()) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no variables.";
    throw err;
  }
  return new MaximumObj(vars.toStdVector());
}

VarObj *
XmlParser::parseVar(QDomElement &node, SimulationObj *sim) {
  if ("var" == node.tagName()) { return parseVarDef(node, sim); }
  if ("ref" == node.tagName()) { return parseVarRef(node, sim); }

  ParserError err;
  err << "ParserError @"<< node.lineNumber()
      << ": Unexpected element " << node.tagName().toStdString()
      << " expected variable definition or reference.";
  throw err;
}

VarObj *
XmlParser::parseVarDef(QDomElement &node, SimulationObj *sim) {
  if (!node.hasAttribute("type")) {
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
  VarObj *var = _factories[node.attribute("type")](node, sim, this);
  if (node.hasAttribute("id")) {
    sim->addVar(node.attribute("id").toStdString(), var);
  }
  return var;
}

VarObj *
XmlParser::parseVarRef(QDomElement &node, SimulationObj *sim) {
  QString name = node.text();
  if (0 == name.size()) {
    ParserError err;
    err << "ParserError @" << node.lineNumber()
        << ": " << node.tagName().toStdString() << " has no content. Expected variable name.";
    throw err;
  }
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
XmlParser::parseParams(QDomElement &node, SimulationObj *sim) {
  QHash<QString, double> params;
  for (QDomElement p=node.firstChildElement("param"); !p.isNull(); p=p.nextSiblingElement("param")) {
    if (! p.hasAttribute("name")) {
      ParserError err;
      err << "ParserError @" << p.lineNumber()
          << ": " << node.tagName().toStdString() << " has no 'name' attribute.";
      throw err;
    }
    bool ok=true; double value = p.text().toDouble(&ok);
    if (! ok) {
      ParserError err;
      err << "ParserError @" << node.lineNumber()
          << ": " << p.tagName().toStdString() << "'s value is not a floating point number.";
      throw err;
    }
    params.insert(p.attribute("name"), value);
  }
  return params;
}

QVector<VarObj *> XmlParser::parseVars(QDomElement &node, SimulationObj *sim) {
  QVector<VarObj *> vars;
  for (QDomElement p=node.firstChildElement(); !p.isNull(); p=p.nextSiblingElement()) {
    vars.push_back(parseVar(p, sim));
  }
  return vars;
}

void XmlParser::parseOutput(QDomElement &node, SimulationObj *sim) {
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
  for (QDomElement p=node.firstChildElement("pdf"); !p.isNull(); p=p.nextSiblingElement("pdf")) {
    QDomElement cp = p.firstChildElement();
    VarObj *var = parseVar(cp, sim);
    sim->addOutputVar(var);
  }
}
