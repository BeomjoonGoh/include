#ifndef LORENTZ_H
#define LORENTZ_H

#include "assert.h"
#include "maths.h"

class Lorentz
{ // Lorentz distribution:
  //  L(\o) = P\frac{\g}{(\o-x_0)^2 + \g^2}
  public:
    double  x0;                                         // Centre of Lorentz distribution
    double  g;                                          // Width of Lorentz distribution, \gamma
    double  P;                                          // Height of Lorentz distribution without pi factor
    bool    exist;                                      // Whether the Lorentz distribution exist or not

  public:
    Lorentz() : x0{0}, g{1}, P{0}, exist{false} { }
    Lorentz(double x0_, double g_, double P_) : x0{x0_}, g{g_}, P{P_}, exist{true} { }
    ~Lorentz() { }

    double operator()(double w);

    void Set(double z0, double e, double a, double p, double q, double r);
    void SetFalse();
    double convLL(const Lorentz& l, double w) const;
    double convLA(double w0, double w1, double A0, double A1, double w) const;
};

inline double Lorentz::operator()(double w)
{
  return P*g/(Maths::sqr(w-x0)+Maths::sqr(g));
}

inline void Lorentz::Set(double z0, double e, double a, double p, double q, double r)
{ // G(\o) ~ \frac{1}{A(\o-z_0)^2 + 2B(\o-z_0) + C} = L(\o) using Taylor expansion
  double A = (Maths::sqr(1-p)+Maths::sqr(q))/a-2*e*q*(r/a)/a+Maths::sqr(e*r/a)/a;
  double B = e*q/a-Maths::sqr(e)/a*(r/a)/2;
  double C = Maths::sqr(e)/a;
  double g2 = C/A-Maths::sqr(B/A);

  x0 = -B/A;
  g = (g2>0)? std::sqrt(g2) : std::sqrt(Maths::abs(C/A));
  if ( g == 0 ) {
    SetFalse();
  } else {
    exist = true;
    x0 += z0;
    P = 1/(A*g);
  }
}

inline void Lorentz::SetFalse()
{
  exist = false;
  P = 0;
}

inline double Lorentz::convLL(const Lorentz& l, double w) const
{ // \int_{-\infty}^{\infty} d\x L(\x+\o) l(\x)
  return P*l.P*Maths::pi*(g+l.g)/(Maths::sqr(g+l.g)+Maths::sqr(x0-l.x0-w));
}


inline double Lorentz::convLA(double w0, double w1, double A0, double A1, double w) const
{ // Analytically treats
  // \int_{\o_0}^{\o_1} d\x \frac{P\g}{(\x+\o-x_0)^2 + \g^2} AF(\x)
  // where, AF(\o) = A_b(+\o)f(-\o), A_b(+\o)f(+\o), or any function
    
  if (!exist)
    return 0.0;
  if (Maths::abs(w1-w0)*100 < g)
    return P*g*0.5*(A0+A1)*(w1-w0)/(Maths::sqr(0.5*(w0+w1)+w-x0)+Maths::sqr(g));

  double c0 = w0 + w - x0;
  double c1 = w1 + w - x0;
  double dA = (A1-A0)/(w1-w0);
  double AdAc = A0 - dA*c0;

  if ( Maths::abs(c0) > 100*g && Maths::abs(c1) > 100*g ) {
    double R = P*g*( AdAc*(1/c0-1/c1) + dA*std::log(Maths::abs(c1/c0)) + 0.5*dA*(Maths::sqr(g/c1)-Maths::sqr(g/c0)) );
    if (c0*c1 > 0.0)     return R;
    if (c1-c0 > 199.9*g) return R + P*AdAc*Maths::pi;
  }

  double a0 = c0/g;
  double a1 = c1/g;

  double R;
  if (Maths::abs(g) < 1e-30)
    R = P*g*( AdAc*(1/c0-1/c1) + dA*std::log(Maths::abs(c1/c0)) );
  else
    R = P*( AdAc*(std::atan(a1)-std::atan(a0)) + 0.5*g*dA*std::log((Maths::sqr(a1)+1)/(Maths::sqr(a0)+1)) );

  assert(!std::isnan(R) && !std::isinf(R),
         "R is nan or inf: "<<R<<" "<< "w0,w1="<<w0<<" "<<w1<<" A0,A1="<<A0<<" "<<A1<<" w,x0="<<w<<" "<<x0<<std::endl
         <<" from: 1+a^2="<<(1+Maths::sqr(a0))<<","<<(1+Maths::sqr(a1))<<" a0="<<a0<<" a1="<<a1<<" g="<<g<<" c0="<<c0<<" c1="<<c1
         <<" std::atan="<<std::atan(a1)-std::atan(a0)<<" AdAc="<<AdAc<<" log"<<std::log((1+Maths::sqr(a1))/(1+Maths::sqr(a0))));
  return R;
}

#endif /* end of include guard: LORENTZ_H */
