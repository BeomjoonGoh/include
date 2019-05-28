#ifndef DILOG_H
#define DILOG_H

#include "maths.h"
#include "complex.h"

class Dilogarithm
{ // Dilogarithm function
  //  Li_2 (z) = \sum_{n=1}^{\infty} \frac{z^n}{n^2}, for |z| < 1.
  // This is a special case of polylogarithm function with m = 2:
  //  Li_m (z) = \sum_{n=1}^{\infty} \frac{z^n}{n^m}, for |z| < 1.
  // Note that for m = 1,
  //  Li_1 (z) = -\log(z-1) = \sum_{n=1}^{\infty} \frac{z^n}{n}, for |z| < 1.
  //
  // B[n] = [Even Bernoulli numbers/(2n+1)]! made by: Table[BernoulliB[2n]/(2n+1)!, {n,1,19}]

  public:
    static double dilog(const double x);
    static compdb dilog(const compdb &z);
  private:
    static const unsigned N = 20;
    static const double C[N];
    static const double B[N];
};

double Dilogarithm::dilog(const double x)
{
  quitif(x > 1, "Li_2(x) is not defined at x=" << x);
  const double pi2 = Maths::sqr(Maths::pi);
  const double pi6 = pi2/6;
  double Y,S,A;
  if (x == 1) {
    return pi6;
  } else if (x == -1) {
    return -pi2/12;
  } else {
    double T = -x;
    if (T <= -2) {
      Y = -1/(1+T);
      S = 1;
      A = -pi2/3+0.5*(Maths::sqr(std::log(-T))-Maths::sqr(std::log(1+1/T)));
    } else if (T < -1) {
      Y = -1-T;
      S = -1;
      A = std::log(-T);
      A = -pi6+A*(A+std::log(1+1/T));
    } else if (T <= -0.5) {
      Y = -(1+T)/T;
      S = 1;
      A = std::log(-T);
      A = -pi6+A*(-0.5*A+std::log(1+T));
    } else if (T < 0) {
      Y = -T/(1+T);
      S = -1;
      A = 0.5*Maths::sqr(std::log(1+T));
    } else if (T <= 1) {
      Y = T;
      S = 1;
      A = 0;
    } else {
      Y = 1/T;
      S = -1;
      A = pi6+0.5*Maths::sqr(std::log(T));
    }
    double H    = Y+Y-1;
    double ALFA = H+H;
    double B2   = 0;
    double B1   = 0;
    double B0;
    for (int i = N-1; i >= 0; i--){
      B0 = C[i] + ALFA*B1-B2;
      B2 = B1;
      B1 = B0;
    }
    return -(S*(B0-H*B2)+A);
  }
}

compdb Dilogarithm::dilog(const compdb& z)
{
  const double pi2 = Maths::sqr(Maths::pi);
  const double rz = real(z);
  const double iz = imag(z);
  const double az = abs(z);

  // special cases
  if (iz == 0.) {
    if (rz <= 1.)
      return compdb{dilog(rz), 0.};
    if (rz > 1.)
      return compdb{dilog(rz), -Maths::pi*std::log(rz)};
  } else if (az < Maths::absEpsilon) {
    return z;
  }

  compdb cy, cx;
  int jsgn, ipi12;
  // transformation to |z|<1, Re(z)<=0.5
  if (rz <= 0.5) {
    if (az > 1.) {
      cy = -0.5 * Maths::sqr(log(-z));
      cx = -log(1. - 1./z);
      jsgn = -1;
      ipi12 = -2;
    } else { // (az <= 1.)
      cy = 0;
      cx = -log(1. - z);
      jsgn = 1;
      ipi12 = 0;
    }
  } else { // rz > 0.5
    if (az <= std::sqrt(2*rz)) {
      cx = -log(z);
      cy = cx * log(1. - z);
      jsgn = -1;
      ipi12 = 2;
    } else { // (az > sqrt(2*rz))
      cy = -0.5 * Maths::sqr(log(-z));
      cx = -log(1. - 1. / z);
      jsgn = -1;
      ipi12 = -2;
    }
  }

  // the dilogarithm
  const compdb cx2{Maths::sqr(cx)};
  compdb sumC;

  // lowest order terms w/ different powers
  for (unsigned i1 = 2; i1 < N; i1++)
    sumC = cx2*(sumC + B[N+1-i1]);

  sumC = cx + cx2*(B[0] + cx*(B[1] + sumC));

  const compdb result = static_cast<double>(jsgn)*sumC + cy + ipi12*pi2/12.;

  return result;
}

const double Dilogarithm::C[20] = {
   0.42996693560813697,
   0.40975987533077105,
  -0.01858843665014592,
   0.00145751084062268,
  -0.00014304184442340,
   0.00001588415541880,
  -0.00000190784959387,
   0.00000024195180854,
  -0.00000003193341274,
   0.00000000434545063,
  -0.00000000060578480,
   0.00000000008612098,
  -0.00000000001244332,
   0.00000000000182256,
  -0.00000000000027007,
   0.00000000000004042,
  -0.00000000000000610,
   0.00000000000000093,
  -0.00000000000000014,
   0.00000000000000002
};
const double Dilogarithm::B[20] = {
  -1./4.,
   1./36.,
  -1./36.e2,
   1./21168.e1,
  -1./108864.e2,
   1./52690176.e1,
  -4.0647616451442255268059093862919666745470571274397078e-11,
   8.9216910204564525552179873167527488515142836130490451e-13,
  -1.9939295860721075687236443477937897056306947496538801e-14,
   4.5189800296199181916504765528555932283968190144666184e-16,
  -1.0356517612181247014483411542218656665960912381686505e-17,
   2.3952186210261867457402837430009803816789490019429743e-19,
  -5.5817858743250093362830745056254199055670546676443981e-21,
   1.3091507554183212858123073991865923017498498387833038e-22,
  -3.0874198024267402932422797648664624315955652561327457e-24,
   7.3159756527022034203579056092521485910334010636908750e-26,
  -1.7408456572340007409890551477597025453408414217542713e-27,
   4.1576356446138997196178996207752266734882541595115639e-29,
  -9.9621484882846221031940067024558388498548600173944888e-31,
   2.3940344248961653005211679878937495629342791569329158e-32,
};

#endif /* end of include guard: DILOG_H */
