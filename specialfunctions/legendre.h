#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "assert.h"
#include "maths.h"

class Legendre
{ // For real x: |x| <= 1, integer m, l: 0 <= m <= l
  // 1) Legendre polynomial P_{l}(x):
  //      \frac{d}{dx} \left[(1-x^2) \frac{d}{dx} P_{l}(x)\right] + l(l+1) P_{l}(x) = 0.
  // 2) Associated Legendre polynomial P^{m}_{l}(x):
  //      P^{m}_{l}(x) = (-1)^{m}(1-x^2)^{m/2} \frac{d^m}{dx^m} P_{l}(x).
  // 3) Renormalized associated Legendre polynomial p^{m}_{l}(x):
  //      p^{m}_{l}(x) = \sqrt{\frac{2l+1}{4\pi} \frac{(l-m)!}{(l+m)!}} P^{m}_{l}(x)
  public:
    static double P(const int l, const double x);
    static double P(const int l, const int m, const double x);
    static double p(const int l, const int m, const double x);
  private:
};

inline double Legendre::P(const int l, const double x)
{ // Legendre polynomial P_{l}(x)
  assert(0 <= l && Maths::abs(x) <= 1.0, "l,x="<<l<<","<<x);
  return std::sqrt(4.0*Maths::pi/(2*l+1))*p(l,0,x);
}

inline double Legendre::P(const int l, const int m, const double x)
{ // Associated Legendre polynomial P^{m}_{l}(x)
  // Note: Will overflow for m ~ 80, or sooner if l >> m
  assert(0 <= m && m <= l && Maths::abs(x) <= 1.0, "m,l,x="<<m<<","<<l<<","<<x);
  double prod = 1.0;
  for (int j = l-m+1; j <= l+m; j++) {
    prod *= j;
  }
  return std::sqrt(4.0*Maths::pi*prod/(2*l+1))*p(l,m,x);
}

inline double Legendre::p(const int l, const int m, const double x)
{ // Renormalized associated Legendre polynomial p^{m}_{l}(x)
  assert(0 <= m && m <= l && Maths::abs(x) <= 1.0, "m,l,x="<<m<<","<<l<<","<<x);
  double pmm = 1.0;
  if (m > 0) {
    double omx2 = (1.0-x)*(1.0+x);
    double fact = 1.0;
    for (int i = 1; i <= m; i++) {
      pmm *= omx2*fact/(fact+1.0);
      fact += 2.0;
    }
  }
  pmm = std::sqrt((2*m+1)*pmm/(4.0*Maths::pi));

  if (m & 1)
    pmm = -pmm;

  if (l == m)
    return pmm;
  else {
    double pmmp1 = x*std::sqrt(2.0*m+3.0)*pmm;
    if (l == m+1)
      return pmmp1;
    else {
      double pll = 0;
      double oldfact = std::sqrt(2.0*m+3.0);
      for (int ll = m+2; ll <= l; ll++) {
        int l2 = Maths::sqr(ll);
        double fact = std::sqrt((4.0*l2-1.0)/(l2-Maths::sqr(m)));
        pll = (x*pmmp1-pmm/oldfact)*fact;

        oldfact = fact;
        pmm = pmmp1;
        pmmp1 = pll;
      }
      return pll;
    }
  }
}

#endif /* end of include guard: LEGENDRE_H */
