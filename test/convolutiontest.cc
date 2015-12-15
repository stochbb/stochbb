#include "convolutiontest.hh"
#include "lib/api.hh"
#include "lib/chain.hh"

using namespace stochbb;
using namespace stochbb::UnitTest;

ConvolutionTest::ConvolutionTest()
{
  // pass...
}

void
ConvolutionTest::testAlignment() {
  Var X1 = normal(100., 10.);
  Var X2 = normal(100., 10.);
  Var X3 = normal(200., std::sqrt(200.));
  Density conv = new ConvolutionDensityObj(std::vector<Var> {X1, X2});

  size_t steps = 1000;
  double Tmin = 0, Tmax = 400;
  Eigen::VectorXd pdf1(steps); conv.eval(Tmin, Tmax, pdf1);
  Eigen::VectorXd pdf2(steps); X3.density().eval(Tmin, Tmax, pdf2);
  for (size_t i=0; i<steps; i++) {
    UT_ASSERT_NEAR_EPS(pdf1(i), pdf2(i), 1e-6);
  }
}

TestSuite *
ConvolutionTest::suite() {
  TestSuite *suite = new TestSuite("Convolution");
  suite->addTest(new TestCaller<ConvolutionTest>("alignment", &ConvolutionTest::testAlignment));
  return suite;
}

