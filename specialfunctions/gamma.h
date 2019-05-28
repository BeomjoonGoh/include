#ifndef GAMMA_H
#define GAMMA_H

#include "assert.h"
#include "maths.h"

struct Gauleg18
{ // Abscissas and weights for Gauss-Legendre quadrature.
    static const int ngau = 18;
    static const double GLy[ngau];
    static const double GLw[ngau];
};

class Gamma : private Gauleg18
{
  public:
    // Gamma function related
    static double logG(const double x);
    static double factorial(const int n);
    static double logFactorial(const int n);
    static double beta(const double z, const double w);

    // Incomplete gamma functions
    static double P(const double a, const double x);
    static double Q(const double a, const double x);
    static double invP(const double p, const double a);

  private:
    static double PseriseExpansion(const double a, const double x);
    static double QcontinuedFraction(const double a, const double x);
    static double Gquadrature(double a, double x, bool p);

  private:
    static const int    BigA = 100;
    static const double Fmin;
};

inline double Gamma::logG(const double x)
{ // Returns the value \ln{\Gamma(x)} for x > 0. \Gamma(z) = \int_0^\infty dt t^{z-1}e^{-t}
  assert(x>0, "x="<<x);
  static const double cof[14] = {
     57.1562356658629235,
    -59.5979603554754912,
     14.1360979747417471,
    -0.491913816097620199,
     3.39946499848118887e-5,
     4.65236289270485756e-5,
    -9.83744753048795646e-5,
     1.58088703224912494e-4,
    -2.10264441724104883e-4,
     2.17439618115212643e-4,
    -1.64318106536763890e-4,
     8.44182239838527433e-5,
    -2.61908384015814087e-5,
     3.68991826595316234e-6
  };
  double y = x;
  double tmp  =  x+5.24218750000000000; // 671/128.
  double fac = (x+0.5)*std::log(tmp)-tmp;
  double ser = 0.999999999999997092;
  for (int j = 0; j < 14; j++)
    ser +=  cof[j]/++y;
  return fac+std::log(2.5066282746310005*ser/x);
}

inline double Gamma::factorial(const int n)
{ // Returns the value n! as a floating-point number.
  static const int NTOP = 171; // 171! overflows.
  static double a[NTOP];
  static bool init = true;
  if (init) {
    init = false;
    a[0] = 1.;
    for (int i = 1; i < NTOP; i++)
      a[i] = i*a[i-1];
  }
  assert(n >= 0 && n < NTOP, "factorial out of range (NTOP=171), n="<<n);
  return a[n];
}

inline double Gamma::logFactorial(const int n)
{ // Returns ln(n!).
  static const int NTOP=2000; // can be increased but necessary?
  static double a[NTOP];
  static bool init = true;
  if (init) {
    init = false;
    for (int i = 0; i < NTOP; i++)
      a[i] = logG(i+1.);
  }
  assert0(n >= 0);
  return (n < NTOP) ? a[n] : logG(n+1.0);
}

inline double Gamma::beta(const double z, const double w)
{ // Returns the value of the beta function B(z,w) = \int_0^1 dt t^{z-1} (1-t)^{w-1}.
  return std::exp(logG(z)+logG(w)-logG(z+w));
}

inline double Gamma::P(const double a, const double x)
{ // Returns the incomplete gamma function P(a,x) = 
  assert(x >= 0.0 && a < 0.0, "x,a="<<x<<","<<a);
  if (x == 0.0)
    return 0.0;
  else if (static_cast<int>(a) >= BigA)
    return Gquadrature(a,x,1);
  else if (x < a+1.0)
    return PseriseExpansion(a,x);
  else
    return 1.0-QcontinuedFraction(a,x);
}

inline double Gamma::Q(const double a, const double x)
{ // Returns the incomplete gamma function Q(a,x) = 1 - P(a,x).
  assert(x >= 0.0 && a < 0.0, "x,a="<<x<<","<<a);
  if (x == 0.0)
    return 1.0;
  else if (static_cast<int>(a) >= BigA)
    return Gquadrature(a,x,0);
  else if (x < a+1.0)
    return 1.0-PseriseExpansion(a,x);
  else
    return QcontinuedFraction(a,x);
}

inline double Gamma::PseriseExpansion(const double a, const double x)
{ // Returns the incomplete gamma function P(a,x) evalutated by its series representation.
  double ap = a;
  double sum = 1.0/a;
  double del = 1.0/a;

  for (;;) {
    ++ap;
    del *= x/ap;
    sum += del;
    if (Maths::abs(del) < Maths::abs(sum)*Maths::absEpsilon)
      return sum*std::exp(-x+a*std::log(x)-logG(a));
  }
}

inline double Gamma::QcontinuedFraction(const double a, const double x)
{ // Returns the incomplete gamma function Q(a,x) evalutated by its continued fraction representation. 
  double b = x + 1.0 - a;
  double c = 1.0/Fmin;
  double d = 1.0/b;
  double fac = d;
  double an, del;
  for (int i = 1;; i++) {
    an = -i*(i-a);
    b += 2.0;
    d = an*d+b;
    if (Maths::abs(d) < Fmin)
      d = Fmin;
    c = b + an/c;
    if (Maths::abs(c) < Fmin)
      c = Fmin;
    d = 1.0/d;
    del = d*c;
    fac *= del;
    if (Maths::abs(del-1.0) <= Maths::absEpsilon) break;
  }
  return fac*std::exp(-x+a*std::log(x)-logG(a));
}

inline double Gamma::Gquadrature(double a, double x, bool p)
{ // Incomplete gamma by quadrature. Returns P(a,x) or Q(a,x), when p = 1 or 0.
  // Directly calcultates the integral by Gauss-Legendre quadrature.
  // Used when a is too large to use PseriseExpansion or QcontinuedFraction
  double a1 = a-1.0;
  double lna1 = std::log(a1);
  double sqrta1 = std::sqrt(a1);
  //Set how far to integrate into the tail:
  double xu = (x > a1) ? Maths::max(a1+11.5*sqrta1, x+6.0*sqrta1) : Maths::max(0.0, Maths::min(a1-7.5*sqrta1, x-5.0*sqrta1));
  double sum = 0;
  double t;
  for (int j = 0; j < ngau; j++) { // Gauss-Legendre.
    t = x + (xu-x)*GLy[j];
    sum += GLw[j]*std::exp(-(t-a1)+a1*(std::log(t)-lna1));
  }
  double ans = sum*(xu-x)*std::exp(a1*(lna1-1.)-logG(a));
  //if (p) return ans > 0.0 ? 1.0-ans : -ans;
  //else   return ans < 0.0 ? 1.0+ans :  ans;
  if (p) return x > a1 ? 1.0-ans : -ans;
  else   return x > a1 ? ans : 1.0+ans;
}

inline double Gamma::invP(const double p, const double a)
{ // Returns x such that P(a,x) = p for an argument p between 0 and 1.
  assert(a > 0.0, "a="<<a);
  if (p >= 1.)
    return Maths::max(100.0, a + 100.0*std::sqrt(a));
  if (p <= 0.)
    return 0.0;

  double x, t;
  double lna1, afac;
  double a1 = a - 1.0;
  double gln = logG(a);
  // initial guess
  if (a > 1.0) {
    lna1 = std::log(a1);
    afac = std::exp(a1*(lna1-1.)-gln);
    double pp = (p < 0.5) ? p : 1.0 - p;
    t = std::sqrt(-2.0*std::log(pp));

    x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
    if (p < 0.5)
      x = -x;
    x = Maths::max(1.0e-3, a*pow(1.0-1.0/(9.0*a)-x/(3.0*std::sqrt(a)),3));
  } else {
    t = 1.0 - a*(0.253+a*0.12);
    x = (p < t) ? pow(p/t,1./a) : 1.0-std::log(1.0-(p-t)/(1.0-t));
  }

  double err, u;
  const double EPS = 1.e-8; // Accuracy is the square of EPS.
  for (int j = 0; j < 12; j++) {
    if (x <= 0.0)
      return 0.0;
    err = P(a,x) - p;
    t = (a > 1.0) ? afac*std::exp(-(x-a1)+a1*(std::log(x)-lna1)) : std::exp(-x+a1*std::log(x)-gln);

    // Halley's method.
    u = err/t;
    t = u/(1.0 - 0.5*Maths::min(1.0, u*(a1/x - 1.0)));
    x -= t;

    if (x <= 0.)
      x = 0.5*(x + t); // Halve old value if x tries to go negative.

    if (Maths::abs(t) < EPS*x ) break;
  }
  return x;
}

const double Gauleg18::GLy[Gauleg18::ngau] = {
  0.0021695375159141994,
  0.011413521097787704,
  0.027972308950302116,
  0.051727015600492421,
  0.082502225484340941,
  0.12007019910960293,
  0.16415283300752470,
  0.21442376986779355,
  0.27051082840644336,
  0.33199876341447887,
  0.39843234186401943,
  0.46931971407375483,
  0.54413605556657973,
  0.62232745288031077,
  0.70331500465597174,
  0.78649910768313447,
  0.87126389619061517,
  0.95698180152629142
};
const double Gauleg18::GLw[Gauleg18::ngau] = {
  0.0055657196642445571,
  0.012915947284065419,
  0.020181515297735382,
  0.027298621498568734,
  0.034213810770299537,
  0.040875750923643261,
  0.047235083490265582,
  0.053244713977759692,
  0.058860144245324798,
  0.064039797355015485,
  0.068745323835736408,
  0.072941885005653087,
  0.076598410645870640,
  0.079687828912071670,
  0.082187266704339706,
  0.084078218979661945,
  0.085346685739338721,
  0.085983275670394821
};

const double Gamma::Fmin = Maths::doubleMin/Maths::absEpsilon;


#endif /* end of include guard: GAMMA_H */
