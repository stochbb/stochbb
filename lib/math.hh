#ifndef __SBB_MATH_HH__
#define __SBB_MATH_HH__

#include <cmath>

namespace sbb {

/** Logarithm of the gamma function.
 * @ingroup internal */
double lgamma(double z);
/** The error function.
 * @ingroup internal */
double erfc(double x);
/** The lower-incomplete gamma function.
 * @ingroup internal */
double gamma_li(double s, double z);
/** The upper-incomplete gamma function.
 * @ingroup internal */
double gamma_ui(double s, double z);
/** The beta function.
 * @ingroup internal */
double betai(double a, double b, double x);

}

#endif // __SBB_MATH_HH__
