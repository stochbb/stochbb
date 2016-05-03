#include "convolutiontest.hh"
#include "lib/api.hh"
#include "lib/chain.hh"
#include "lib/distribution.hh"


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
  Var X4 = X1 + X2;
  Density conv = new ConvolutionDensityObj(std::vector<Var> {X1, X2});

  size_t steps = 1000;
  double Tmin = 0, Tmax = 400;
  Eigen::VectorXd pdf1(steps); conv.eval(Tmin, Tmax, pdf1);
  Eigen::VectorXd pdf2(steps); X3.density().eval(Tmin, Tmax, pdf2);
  Eigen::VectorXd pdf3(steps); X4.density().eval(Tmin, Tmax, pdf3);
  for (size_t i=0; i<steps; i++) {
    UT_ASSERT_NEAR_EPS(pdf1(i), pdf2(i), 1e-6);
    UT_ASSERT_NEAR_EPS(pdf3(i), pdf2(i), 1e-6);
  }
}

void
ConvolutionTest::testDeltaReduction() {
  Var X = delta(10);
  Var Y = gamma(10,10);
  Var Z = X + Y;

  UT_ASSERT(Z.density().is<AtomicDensity>());
  // Cast to atomic density
  AtomicDensity Zatom = Z.density().as<AtomicDensity>();
  // Check distribution type
  UT_ASSERT(Zatom.distribution().is<GammaDistribution>());
  // Check parameters
  UT_ASSERT_EQUAL(10., Zatom.parameter(0));
  UT_ASSERT_EQUAL(10., Zatom.parameter(1));
  UT_ASSERT_EQUAL(10., Zatom.parameter(2));
}

void
ConvolutionTest::testNormalReduction() {
  Var X = normal(1,1);
  Var Y = normal(1,1);
  Var Z = X + Y;

  UT_ASSERT(Z.density().is<AtomicDensity>());
  // Cast to atomic density
  AtomicDensity Zatom = Z.density().as<AtomicDensity>();
  // Check distribution type
  UT_ASSERT(Zatom.distribution().is<NormalDistribution>());
  // Check parameters
  UT_ASSERT_EQUAL(2., Zatom.parameter(0));
  UT_ASSERT_EQUAL(std::sqrt(2), Zatom.parameter(1));
}

void
ConvolutionTest::testGammaReduction() {
  Var X = gamma(1,1);
  Var Y = gamma(1,1);
  Var Z = X + Y;

  UT_ASSERT(Z.density().is<AtomicDensity>());
  // Cast to atomic density
  AtomicDensity Zatom = Z.density().as<AtomicDensity>();
  // Check distribution type
  UT_ASSERT(Zatom.distribution().is<GammaDistribution>());
  // Check parameters
  UT_ASSERT_EQUAL(2., Zatom.parameter(0));
  UT_ASSERT_EQUAL(1., Zatom.parameter(1));
  UT_ASSERT_EQUAL(0., Zatom.parameter(2));
}

TestSuite *
ConvolutionTest::suite() {
  TestSuite *suite = new TestSuite("Convolution");
  suite->addTest(new TestCaller<ConvolutionTest>("alignment", &ConvolutionTest::testAlignment));
  suite->addTest(new TestCaller<ConvolutionTest>("delta reduction", &ConvolutionTest::testDeltaReduction));
  suite->addTest(new TestCaller<ConvolutionTest>("normal reduction", &ConvolutionTest::testNormalReduction));
  suite->addTest(new TestCaller<ConvolutionTest>("gamma reduction", &ConvolutionTest::testNormalReduction));
  return suite;
}

