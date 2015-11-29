#include "lib/randomvariable.hh"
#include "lib/operators.hh"
#include "lib/minmax.hh"
#include "lib/xmlparser.hh"
#include "lib/exception.hh"

#include <iostream>
#include <QFile>
#include <QDomDocument>

using namespace sbb;

int main(int argc, char *argv[]) {
  if (2 > argc) {
    std::cerr << "Usage: stochbb FILENAME" << std::endl;
    return -1;
  }

  QFile file(argv[1]);
  if (! file.open(QIODevice::ReadOnly)) {
    std::cerr << "Cannot open file " << argv[1] << std::endl;
    return -1;
  }

  QDomDocument doc;
  QString msg; int row;
  if (! doc.setContent(&file, &msg, &row)) {
    file.close();
    std::cerr << "Cannot parse file " << argv[1] << ": " << msg.toStdString()
              << " @" << argv[1] << ":" << row << std::endl;
    return -1;
  }

  XmlParser parser;
  Simulation sim;
  QDomElement root = doc.documentElement();
  try { sim = parser.parse(root); }
  catch (Error &err) {
    std::cerr << err.what() << std::endl;
    return -1;
  }

  size_t N = sim.steps(), M = sim.numOutputVars();
  Eigen::MatrixXd out; sim.run(out);
  for (size_t i=0; i<N; i++) {
    std::cout << out(i,0);
    for (size_t j=0; j<M; j++) {
      std::cout << '\t' << out(i,j+1);
    }
    std::cout << std::endl;
  }

  return 0;
}
