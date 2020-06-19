#ifndef DERULE_H
#define DERULE_H

#include <iostream>
#include "assert.h"
#include "maths.h"

namespace Maths {
  template <class functor> double integral(double a, double b, functor &f, double tol = absEpsilon);
  template <class functor> double integral(double a, double b, functor &f, int m);
}

template <class functor>
class DErule
{ //  I = \int_a^b dx f(x)
  // Assumptions: 1) Interior is smooth enough.
  //              2) f(a) or f(b) can be at least an integrable singularity.
  //
  // Quadrature by variable transformation using DE (Double Exponential) rule:
  //  I = \int_{-\infty}^{\infty} dt f(x(t)) x'(t)
  // where,
  //  x(t)  = \frac{1}{2}(b+a) + \frac{1}{2}(b-a) \tanh(\sinh(t)),
  //  x'(t) = \frac{dx}{dt} = \frac{1}{2}(b-a) \sech^{2}(\sinh(t))\cosh(t).
  //
  // The quadrature is calculated as (with \d_j = b - x(t_j), t_{\pmN} = \pm hmax),
  //  I ~ h \sum_{j=-N}^{N} f(x_j) x'(t_j)
  //    = h \left\{ f((a+b)/2) x'(t_0) + \sum_{j=1}^{N} [ f(a+\d_j) + f(b-\d_j) ] x'(t_j) \right\}
  // Note:
  //  1)  x'(t) ~ \exp(\exp|t|) as |t| \to \infty, hence the name double exponential
  //  2)  This is trapezoidal rule with no 1/2 factors at ends; they become irrelevant due to DE transformation.
  //  3)  It's written in symmetric form to handle overflow from \sech^{2}(\sinh t). This is done by defining
  //      q = \exp(-2\sinh t) (q underflows for positive t, and the symmetric form handles negative t).
  //  4)  The function/functor f(x) requires two double arguments: (double x, double delta);
  //      If the singularities are no worse than logarithm, no delta adjusted functor is needed (treat delta as dummy).
  //  5)  Use hmax = 3.7 for mild singularities, higher value for stronger singularities.

  private:
    functor &f;         // double f(double, double)
    const double a, b;  // Limits of integration
    const double hmax;  // Limits of integration for t

    int n;              // Current level of refinement
    double sum;         // Current value of integral

  public:
    DErule(functor &f_, double a_, double b_, double H = 3.7) : f(f_), a{a_}, b{b_}, hmax{H}, n{0}, sum{0.0} { }
    ~DErule() { }

    int    level() { return n; }
    double value() { return sum; }
    void   next();
};

template <class functor>
inline void DErule<functor>::next()
{
  n++;
  if (n == 1) {
    sum = hmax*0.5*(b-a)*f(0.5*(b+a), 0.5*(b-a));
  } else {
    int it = 1<<(n-2);
    double twoh = hmax/it;
    double t = 0.5*twoh;
    double ssum = 0.0;
    for (int j = 0; j < it; j++) {
      double q = exp(-2.0*sinh(t));
      double del = (b-a)*q/(1.0+q);
      double fac = q/Maths::sqr(1.0+q)*cosh(t);
      ssum += fac*(f(a+del, del) + f(b-del, del));
      t += twoh;
    }
    sum = 0.5*sum + (b-a)*twoh*ssum;
  }
}

template <class functor>
inline double Maths::integral(double a, double b, functor &f, double tol)
{ // returns \int_a^b dx f(x) using DE rule, and recursively halving abscissa with next() function until convergence
  if (abs(a - b) < tol)
    return 0.0;

  const int nmax = 18; // Maximum step.
  const int nmin = 3;  // Warm-up step.

  DErule<functor> I(f,a,b,3.7);

  while (I.level() <= nmin)
    I.next();
  double old = I.value();
  while (I.level() <= nmax) {
    I.next();
    if (abs(I.value()-old) < tol*abs(old) || (I.value() == 0.0 && old == 0.0)) {
      break;
    }
    old = I.value();
  }

  quitif(I.level() > nmax, RED("Not converged!")<<"(nmax=18), I.level()="<<I.level());
  return I.value();
}

template <class functor>
inline double Maths::integral(double a, double b, functor &f, int m)
{ // returns \int_a^b dx f(x) using DE rule, and recursively halving abscissa with next() function for m times
  DErule<functor> I(f,a,b);
  while (I.level() <= m)
    I.next();
  return I.value();
}

#endif /* end of include guard: DERULE_H */
