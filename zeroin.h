#ifndef ZEROIN_H
#define ZEROIN_H

/**********************************************************************************************************************
 *
 *                        Numerical Math Package
 *
 *                          Brent's root finder
 *             obtains a zero of a function of one variable
 *
 * Synopsis
 *      double zeroin(ax,bx,f,tol=EPSILON)
 *      const double ax                 The root is to be sought within
 *      const double bx                 the interval [ax,bx]
 *      UnivariateFunctor& f            The function under consideration
 *      const double tol                Acceptable tolerance for the root position. It is an optional parameter with
 *                                      default value dbl_epsilon
 *
 *      Zeroin returns an approximate location for the root with accuracy 4*dbl_epsilon*abs(x) + tol
 *
 * Algorithm
 *      G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical computations.
 *      M., Mir, 1980, p.180 of the Russian edition
 *
 * The function makes use of a bissection procedure combined with a linear or quadratic inverse interpolation.
 * At each step the code operates three abscissae - a, b, and c:
 *      b - the last and the best approximation to the root
 *      a - the last but one approximation
 *      c - the last but one or even an earlier approximation such that
 *              1) |f(b)| <= |f(c)|
 *              2) f(b) and f(c) have opposite signs, i.e. b and c encompass the root
 *
 * Given these abscissae, the code computes two new approximations, one by the bissection procedure and the other one
 * from interpolation (if a,b, and c are all different the quadratic interpolation is used, linear otherwise).  If the
 * approximation obtained by the interpolation looks reasonable (i.e. falls within the current interval [b,c], not too
 * close to the end points of the interval), the point is accepted as a new approximation to the root. Otherwise, the
 * result of the bissection is used. Therefore, the range of uncertainty is guaranteed to tighten at least by a factor
 * of 1.6
 *
 *********************************************************************************************************************/

#include "assert.h"

template <class functor>
inline double zeroin(const double ax, const double bx, functor& f, const double tol)
{
  double dbl_epsilon = 2.22045e-16;
  assert(tol > 0, "tol="<<tol);
  assert(bx > ax, "Left end point of the interval should be strictly less than the right one: bx,ax="<<bx<<","<<ax);

  double b = bx;              // the last and the best approx to the root
  double fb = f(b);
  double a = ax;              // the last but one approximation
  double fa = f(a);
  double c = a;               // the last but one or even an earlier approx
  double fc = fa;             // see the condition above

  for(;;) {
    const double prev_step = b-a;           // Step from the previous iteration

    if ( std::fabs(fc) < std::fabs(fb) ) {    // Swap data so that b would be thebest approximation found so far
      a = b;  b = c;  c = a;
      fa=fb;  fb=fc;  fc=fa;
    }
    // Estimate the effective tolerance
    const double tol_act = 2*dbl_epsilon*std::fabs(b) + tol/2;
    double new_step = (c-b)/2;                              // Bissection step for this iteration

    if ( std::fabs(new_step) <= tol_act || fb == 0 )
      // Acceptable approximation is found
      return b;

    // Figuring out if the interpolation can be tried
    if ( std::fabs(prev_step) >= tol_act && std::fabs(fa) > std::fabs(fb) ) {
      // If prev_step was large enough and was in true direction, Interpolation may be tried

      double p;           // Interpolation step is calculated in the form p/q;
      double q;           // division operations is delayed until the last moment
      const double cb = c-b;

      if ( a==c ) {       // If we've got only two distinct points linear interpolation can only be applied
        register const double t1 = fb/fa;
        p = cb*t1;
        q = 1.0 - t1;
      } else {            // Quadratic inverse interpolation
        register const double t1=fb/fc, t2=fb/fa;
        q = fa/fc;
        p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
        q = (q-1.0) * (t1-1.0) * (t2-1.0);
      }

      if ( p > 0 )        // Formulas above computed new_step = p/q with the wrong sign (on purpose).
        q = -q;
      else                // Correct this, but in such a way so that p would be positive
        p = -p;

      if ( 2*p < (1.5*cb*q-std::fabs(tol_act*q)) && 2*p < std::fabs(prev_step*q) )
        // If b+p/q falls in [b,c]  and isn't too large it is accepted
        new_step = p/q;

      // If p/q is too large then the bissection procedure can reduce [b,c] to a larger extent
    }

    if ( std::fabs(new_step) < tol_act ) 
      // Adjust the step to be not less than the tolerance
      new_step =  new_step > 0 ?  tol_act : -tol_act;

    a = b;  fa = fb;                        // Save the previous approximation
    b += new_step;  fb = f(b);              // Do step to a new approximation

    if ( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
      // Adjust c for it to have the sign opposite to that of b
      c = a;  fc = fa;
    }
  }
}

#endif /* end of include guard: ZEROIN_H */
