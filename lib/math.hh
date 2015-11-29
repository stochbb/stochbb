#ifndef __SBB_MATH_HH__
#define __SBB_MATH_HH__

#include <cmath>

namespace sbb {

double lgamma(double z);
double erfc(double x);
double gamma_li(double s, double z);
double gamma_ui(double s, double z);
double betai(double a, double b, double x);

}

#endif // __SBB_MATH_HH__
