#ifndef PLOT_HH
#define PLOT_HH

#include <Eigen/Eigen>
#include <vector>


void do_plot(int argc, char *argv[], const Eigen::MatrixXd &out, const std::vector<std::string> &names);
void output_csv(Eigen::MatrixXd &out, std::ostream &stream);

#endif // PLOT_HH
