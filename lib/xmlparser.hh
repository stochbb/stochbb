#ifndef __SBB_XMLPARSER_HH__
#define __SBB_XMLPARSER_HH__

#include "api.hh"

#include <QDomElement>
#include <QVector>
#include <QHash>

namespace sbb {

class XmlParser
{
public:
  XmlParser();
  Simulation parse(const QString &filename);
  Simulation parse(QDomElement &doc);

protected:
  VarObj *parseVar(QDomElement &node, SimulationObj *sim);
  VarObj *parseVarDef(QDomElement &node, SimulationObj *sim);
  VarObj *parseVarRef(QDomElement &node, SimulationObj *sim);
  void parseOutput(QDomElement &node, SimulationObj *sim);

  static VarObj *parseDelta(QDomElement &node, SimulationObj *sim, XmlParser *parser);
  static VarObj *parseUnif(QDomElement &node, SimulationObj *sim, XmlParser *parser);
  static VarObj *parseNorm(QDomElement &node, SimulationObj *sim, XmlParser *parser);
  static VarObj *parseGamma(QDomElement &node, SimulationObj *sim, XmlParser *parser);
  static VarObj *parseChain(QDomElement &node, SimulationObj *sim, XmlParser *parser);
  static VarObj *parseMaximum(QDomElement &node, SimulationObj *sim, XmlParser *parser);
  static VarObj *parseMinimum(QDomElement &node, SimulationObj *sim, XmlParser *parser);

  QHash<QString, double> parseParams(QDomElement &node, SimulationObj *sim);
  QVector<VarObj *> parseVars(QDomElement &node, SimulationObj *sim);

protected:
  QHash<QString, VarObj *(*)(QDomElement &, SimulationObj *, XmlParser *)> _factories;
};

}

#endif // __SBB_XMLPARSER_HH__
