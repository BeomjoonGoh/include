#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "assert.h"
#include "maths.h"

class Ellip
{ // Elliptic integrals
  // Legendre elliptic integrals of the first, second, third kind:
  //    F(\phi,k)     = \int_0^\phi \frac{d\theta}{sqrt{1-k^2\std::sin^2\theta}}
  //    E(\phi,k)     = \int_0^\phi d\theta sqrt{1-k^2\std::sin^2\theta},
  //    \Pi(\phi,n,k) = \int_0^\phi \frac{d\theta}{(1+n\std::sin^2\theta)sqrt{1-k^2\std::sin^2\theta}}
  // (valid for 0 <= \phi <= \pi/2, 0 <= k\std::sin\phi <= 1. Note the sign convention on n.)
  //
  // Complete elliptic integrals of the first, second, third kind (\phi = \pi/2)
  //    K(k)     = \int_0^{\pi/2} \frac{d\theta}{sqrt{1-k^2\std::sin^2\theta}}
  //    E(k)     = \int_0^{\pi/2} d\theta sqrt{1-k^2\std::sin^2\theta}
  //    \Pi(n,k) = \int_0^{\pi/2} \frac{d\theta}{(1+n\std::sin^2\theta)sqrt{1-k^2\std::sin^2\theta}}
  // (valid for 0 <= k <= 1. Note the sign convention on n.)
  //
  // Carlson's symmetric form of elliptic integrals of the first, second, third kind and the degenerate form
  //    R_F(x,y,z)   = \frac{1}{2}\int^\infty_0 \frac{dt}{\sqrt{(t+x)(t+y)(t+z)}}
  //    R_D(x,y,z)   = \frac{3}{2}\int^\infty_0 \frac{dt}{(t+z)\sqrt{(t+x)(t+y)(t+z)}}  = R_J(x,y,z,z)
  //    R_J(x,y,z,p) = \frac{3}{2}\int^\infty_0 \frac{dt}{(t+p)\sqrt{(t+x)(t+y)(t+z)}}
  //    R_C(x,y)     = \frac{1}{2}\int^\infty_0 \frac{dt}{(t+y)\sqrt{(t+x)}}            = R_F(x,y,y)

  public:
    static double legendre1(const double phi, const double k);
    static double legendre2(const double phi, const double k);
    static double legendre3(const double phi, const double n, const double k);

    static double complete1(const double k);
    static double complete2(const double k);
    static double complete3(const double n, const double k);
  
    static double Rf(const double x, const double y, const double z);
    static double Rd(const double x, const double y, const double z);
    static double Rj(const double x, const double y, const double z, const double p);
    static double Rc(const double x, const double y);
};

double Ellip::legendre1(const double phi, const double k)
{ // F(\phi,k) evaluated using R_F.
  double s = std::sin(phi);
  double c = Maths::sqr(std::cos(phi));
  double q = (1.0-s*k)*(1.0+s*k);
  return s*Rf(c, q, 1.0);
}

double Ellip::legendre2(const double phi, const double k)
{ // E(\phi,k) evaluated using R_D, R_F.
  double s = std::sin(phi);
  double c = Maths::sqr(std::cos(phi));
  double q = (1.0-s*k)*(1.0+s*k);
  return s*(Rf(c, q, 1.0) - (Maths::sqr(s*k)/3.0)*Rd(c, q, 1.0));
}

double Ellip::legendre3(const double phi, const double n, const double k)
{ // \Pi(\phi,n,k) evaluated using R_J, R_F.
  double s = std::sin(phi);
  double c = Maths::sqr(std::cos(phi));
  double q = (1.0-s*k)*(1.0+s*k);
  double ns = n*Maths::sqr(s);
  return s*(Rf(c, q, 1.0) - (ns/3.0)*Rj(c, q, 1.0, 1.0+ns));
}

double Ellip::complete1(const double k)
{ // K(k) evaluated using R_F.
  double k2 = Maths::sqr(k);
  return Rf(0.0, 1.0-k2, 1.0);
}

double Ellip::complete2(const double k)
{ // E(k) evaluated using R_D, R_F.
  double k2 = Maths::sqr(k);
  return Rf(0.0, 1.0-k2, 1.0) - (k2/3.0)*Rd(0.0, 1.0-k2, 1.0);
}

double Ellip::complete3(const double n, const double k)
{ // \Pi(n,k) evaluated using R_J, R_F.
  double k2 = Maths::sqr(k);
  return Rf(0.0, 1.0-k2, 1.0) - (n/3.0)*Rj(0.0, 1.0-k2, 1.0, 1.0+n);
}

double Ellip::Rf(const double x_, const double y_, const double z_)
{ // R_F(x, y, z). x, y, and z must be non-negative, and at most one can be zero
  static const double Tol = 0.0025; // for double precision. for single precision, use 0.08
  static const double C[4] = {
    -0.10000000000000000, // -1/10
     0.04166666666666667, //  1/24
    -0.06818181818181818, // -3/44
     0.07142857142857142  //  1/14
  };
  assert(Maths::min(Maths::min(x_,y_),z_) >= 0.0, "x,y,z="<<x_<<","<<y_<<","<<z_<<" must be non-negative.");
  assert(Maths::doubleMin*5.0 <= Maths::min(Maths::min(x_+y_,x_+z_),y_+z_)
         && Maths::max(Maths::max(x_,y_),z_) <= Maths::doubleMax/5.0,
         "Arguments must be in range of double.");

  double x = x_, y = y_, z = z_;
  double mu, X, Y, Z;
  do {
    double sx = std::sqrt(x), sy = std::sqrt(y), sz = std::sqrt(z);
    double l = sx*(sy+sz) + sy*sz;
    x = (x+l)/4.0;
    y = (y+l)/4.0;
    z = (z+l)/4.0;

    mu = (x+y+z)/3.0;
    X = 1.0-x/mu;
    Y = 1.0-y/mu;
    Z = 1.0-z/mu;
  } while (Maths::max(Maths::max(Maths::abs(X),Maths::abs(Y)),Maths::abs(Z)) > Tol);
  double Ea = X*Y - Maths::sqr(Z);
  double Eb = X*Y*Z;
  return (1.0 + (C[0] + C[1]*Ea + C[2]*Eb)*Ea + C[3]*Eb)/std::sqrt(mu);
}

double Ellip::Rd(const double x_, const double y_, const double z_)
{ // R_D(x, y, z). x, y must be non-negative, and at most one can be zero. z must be positive.
  static const double Tol = 0.0015; // for double precision. for single precision, use 0.05
  static const double C[6] = {
    -0.21428571428571427, // -3/14
     0.16666666666666667, //  1/6
    -0.40909090909090909, // -9/22
     0.11538461538461539, //  3/26
     0.10227272727272727, //  9/88
    -0.17307692307692307  // -9/52
  };
  assert(Maths::min(x_,y_) >= 0.0, "x,y="<<x_<<","<<y_<<" must be non-negative.");
  assert(2.0*std::pow(Maths::doubleMax,-2.0/3.0) <= Maths::min(x_+y_,z_)
         && Maths::max(Maths::max(x_,y_),z_) <= 0.1*Tol*std::pow(Maths::doubleMin,-2.0/3.0),
         "Arguments must be in range of double.");

  double x = x_, y = y_, z = z_;
  double mu, X, Y, Z;
  double sum = 0.0;
  double fac = 1.0;
  do {
    double sx = std::sqrt(x), sy = std::sqrt(y), sz = std::sqrt(z);
    double l = sx*(sy+sz) + sy*sz;
    sum += fac/(sz*(z+l));
    fac /= 4.0;
    x = (x+l)/4.0;
    y = (y+l)/4.0;
    z = (z+l)/4.0;

    mu = (x+y+3.0*z)/5.0;
    X = 1.0-x/mu;
    Y = 1.0-y/mu;
    Z = 1.0-z/mu;
  } while (Maths::max(Maths::max(Maths::abs(X),Maths::abs(Y)),Maths::abs(Z)) > Tol);
  double Ea = X*Y;
  double Eb = Maths::sqr(Z);
  double Ec = Ea-Eb;
  double Ed = Ea-6.0*Eb;
  double Ee = Ed+2.0*Ec;
  return 3.0*sum + fac*(1.0 + (C[0]+C[4]*Ed+C[5]*Z*Ee)*Ed + Z*(C[1]*Ee+Z*(C[2]*Ec+Z*C[3]*Ea)))/(mu*std::sqrt(mu));
}

double Ellip::Rj(const double x_, const double y_, const double z_, const double p_)
{ // R_J(x, y, z, p). x, y and z must be non-negative, and at most one can be zero. p must be nonzero.
  // If p < 0, the Cauchy principal value is returned.
  static const double Tol = 0.0015; // for double precision. for single precision, use 0.05
  static const double C[8] = {
    -0.21428571428571427, // -3/14
     0.33333333333333333, //  1/3
    -0.13636363636363636, // -3/22
     0.11538461538461539, //  3/26
     0.10227272727272727, //  9/88
    -0.17307692307692307, // -9/52
     0.16666666666666667, //  1/6
    -0.27272727272727272  // -3/11
  };
  assert(Maths::min(Maths::min(x_,y_),z_) >= 0.0, "x,y,z="<<x_<<","<<y_<<","<<z_<<" must be non-negative.");
  assert(std::pow(Maths::doubleMin*5.0, 1.0/3.0) <= Maths::min(Maths::min(x_+y_,x_+z_), Maths::min(y_+z_,Maths::abs(p_)))
         && Maths::max(Maths::max(x_,y_),Maths::max(z_,Maths::abs(p_))) <= 0.3*std::pow(Maths::doubleMax/5.0, 1.0/3.0),
         "Arguments must be in range of double.");

  double x, y, z, p;
  double a, b, Rcx;
  if (p_ > 0.0) {
    x = x_;
    y = y_;
    z = z_;
    p = p_;
    a = b = Rcx = 0.0;
  } else {
    x = Maths::min(Maths::min(x_,y_),z_);
    z = Maths::max(Maths::max(x_,y_),z_);
    y = x_+y_+z_-x-z;
    a = 1.0/(y-p_);
    b = a*(z-y)*(y-x);
    p = y+b;
    Rcx = Rc(x*z/y,p_*p/y);
  }

  double mu, X, Y, Z, P;
  double sum = 0.0;
  double fac = 1.0;
  do {
    double sx = std::sqrt(x), sy = std::sqrt(y), sz = std::sqrt(z);
    double l = sx*(sy+sz) + sy*sz;

    double alpha = Maths::sqr(p*(sx+sy+sz)+sx*sy*sz);
    double beta = p*Maths::sqr(p+l);
    sum += fac*Rc(alpha,beta);
    fac /= 4.0;
    x = (x+l)/4.0;
    y = (y+l)/4.0;
    z = (z+l)/4.0;
    p = (p+l)/4.0;

    mu = (x+y+z+2.0*p)/5.0;
    X = 1.0-x/mu;
    Y = 1.0-y/mu;
    Z = 1.0-z/mu;
    P = 1.0-p/mu;
  } while (Maths::max(Maths::max(Maths::abs(X),Maths::abs(Y)), Maths::max(Maths::abs(Z),Maths::abs(P))) > Tol);
  double Ea = X*(Y+Z)+Y*Z;
  double Eb = X*Y*Z;
  double Ec = Maths::sqr(P);
  double Ed = Ea-3.0*Ec;
  double Ee = Eb+2.0*P*(Ea-Ec);

  double R = 3.0*sum + fac*(1.0 + Ea*P*(C[1] + P*C[2]) + Eb*(C[6] + P*(C[7] + P*C[3]))
                                - Ec*P*C[1] + Ed*(C[0] + C[4]*Ed + C[5]*Ee))/(mu*std::sqrt(mu));
  if (p <= 0.0) 
    R = a*(b*R + 3.0*(Rcx-Rf(x,y,z)));
  return R;
}

double Ellip::Rc(const double x_, const double y_)
{ // R_C(x, y). x must be non-negative and y must be nonzero. If y < 0, the Cauchy principal value is returned.
  static const double Tol = 0.0012; // for double precision. for single precision, use 0.04
  static const double C[4] = {
     0.30000000000000000, //  3/10
     0.14285714285714285, //  1/7
     0.37500000000000000, //  3/8
     0.40909090909090909  //  9/22
  };
  assert(x_ >= 0.0 && y_ != 0.0, "x_,y_="<<x_<<","<<y_);
  assert(x_+Maths::abs(y_) >= Maths::doubleMin*5.0 && x_+Maths::abs(y_) <= Maths::doubleMax/5.0
         && !(y_ < -std::sqrt(Maths::doubleMin) && x_ > 0.0 && x_ < Maths::sqr(Maths::doubleMin*Maths::doubleMax/5.0)),
         "Arguments must be in range of double.");

  double x, y, w;
  if (y_ > 0.0) {
    x = x_;
    y = y_;
    w = 1.0;
  } else {
    x = x_-y_;
    y = -y_;
    w = std::sqrt(x_)/std::sqrt(x);
  }

  double mu, s;
  do {
    double l = 2.0*std::sqrt(x)*std::sqrt(y)+y;
    x = (x+l)/4.0;
    y = (y+l)/4.0;

    mu = (x+2.0*y)/3.0;
    s = y/mu - 1.0;
  } while (Maths::abs(s) > Tol);

  return w*(1.0+s*s*(C[0]+s*(C[1]+s*(C[2]+s*C[3]))))/std::sqrt(mu);
}


#endif /* end of include guard: ELLIPTIC_H */
