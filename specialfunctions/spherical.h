#ifndef SPHERICAL_H
#define SPHERICAL_H

#include "assert.h"
#include "maths.h"
#include "complex.h"
#include "legendre.h"

class Spherical
{ // Spherical harmonics Y_{lm}(\theta, \phi), -l <= m <= l.
  // they are related to the renormalized associated Legendre polynomials for m >= 0:
  //    Y_{lm}(\theta,\phi) = p^{m}_{l}(cos\theta) e^{im\phi}
  // and 
  //    Y_{l-m}(\theta,\phi) = (-1)^m [Y_{lm}(\theta,\phi)]^{*}
  // so for m < 0 as well.
  public:
    static compdb Y(const int l, const int m, const double theta, const double phi)
    {
      int M = std::abs(m);
      assert(0 <= l && M <= l, "l,M="<<l<<","<<M);

      double plMx = Legendre::p(l,M,std::cos(theta));
      if (m < 0 && M & 1)
        plMx = -plMx;

      double mphi = m*phi;
      return compdb{plMx*std::cos(mphi), plMx*std::sin(mphi)};
    }
};


#endif /* end of include guard: SPHERICAL_H */
