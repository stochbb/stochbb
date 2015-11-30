#include "lib/api.hh"
#include "lib/option_parser.hh"
#include "plot.hh"

#include <iostream>
#include <fstream>

using namespace sbb;

int main(int argc, char *argv[])
{
  /*
   * Assemble command-line options grammar...
   */
  opt::Parser parser;
  opt::RuleInterface &cmd_help
      = parser.Flag("help");
  opt::RuleInterface &cmd_version
      = parser.Flag("version");
  opt::RuleInterface &cmd_run =
      (parser.opt(parser.Flag("log-debug")),
       parser.Value("filename"),
       parser.opt(parser.Flag("plot") | parser.Option("csv")));
  parser.setGrammar( cmd_help | cmd_version | cmd_run );

  // parse command line options
  if (! parser.parse((const char **)argv, argc)) {
    std::cerr << "Error while parsing arguments." << std::endl
              << parser.format_help("stochbb") << std::endl;
    return -1;
  }

  // on "help" command
  if (parser.has_flag("help")) {
    std::cout << parser.format_help("wtcli")  << std::endl;
    return 0;
  }

  // on "version" command
  if (parser.has_flag("version")) {
    std::cout << "stochbb - version 0.0" << std::endl;
    return 0;
  }

  IOLogHandler *handler = 0;
  if (parser.has_flag("log-debug")) {
    handler = new IOLogHandler(std::cerr, LogMessage::DEBUG);
  } else {
    handler = new IOLogHandler(std::cerr, LogMessage::INFO);
  }
  Logger::addHandler(handler);

  // Otherwise run "simulation"
  Simulation sim;
  // get & parse filename
  std::string filename = parser.get_values("filename").front();
  logDebug() << "Load file '" << filename << "'.";
  try { sim = Simulation::fromXml(filename); }
  catch (Error &err) {
    std::cerr << "Error while parsing file " << filename
              << ": " << err.what() << std::endl;
    return -1;
  }

  // run simulation
  size_t N = sim.steps(), M = sim.numOutputVars();
  Eigen::MatrixXd out(N, M+1); sim.run(out);

  // assemble vector of output variable names
  std::vector<std::string> names; names.reserve(M);
  for (size_t i=0; i<M; i++) {
    names.push_back(sim.outputVar(i).name());
  }

  if (parser.has_flag("plot")) {
    // if plot flag is specified -> plot results
    do_plot(argc, argv, out, names);
  } else if (parser.has_option("csv")) {
    std::string filename = parser.get_option("csv").front();
    std::ofstream file;
    file.open(filename.c_str());
    if (! file.is_open()) {
      std::cerr << "Cannot open file '" << filename << "' for output." << std::endl;
      return -1;
    }
    output_csv(out, file);
    file.flush(); file.close();
  } else {
    // otherwise dump as csv
    output_csv(out, std::cout);
  }

  return 0;
}
