#ifndef NONINTDOS_H
#define NONINTDOS_H

#include <iostream>

#include "maths.h"
#include "function.h"
#include "mesh.h"
#include "elliptic.h"
#include "derule.h"

#include <omp.h>

using namespace std;
using namespace Maths;

class Dos
{ // Non interacting density of states of simple cupic lattice for 1, 2, and 3 dimension.
  // E. N. Economou, Green's Functions in Quantum Physics, 3rd Ed., (Springer, 2006), p. 87.
  // T. Morita and T. Horiguchi, J. Math. Phys. 12 981 (1971).
  private:
    const int d;
    const int N;
    const int n;
    const double one;
    double *T;
    F1d<double> D;
  public:
    M1d ep;

  public:
    Dos(int d_, int N_, double t);
    virtual ~Dos() { delete[] T; }

    double& operator[] (const int i) { return D[i]; }
    const double& operator[] (const int i) const { return D[i]; }
    double de(int i) const
    {
      if (i == 0 || i == N-1) {
        return ep.dx(i)*2;
      }
      return ep.dx(i);
    }
    double d1e(int i) const { return ep.d1x(i); }

    int size() const { return N; }
    void setDos();
    double RealG0(double om);
    double ImagG0(double om);
    void Print(const string &outputf, const string &comment = "") { D.Print(ep,outputf,comment); } 

  private:
    double integral0(double a, double b, double e);
    double integral1(double a, double b, double e);
    double integral2(double a, double b, double e);
    double integral3(double a, double b, double e);
};

Dos::Dos(int d_, int N_, double t) : d{d_}, N{N_>>1<<1}, n{3*N}, one{1.0-absEpsilon}, D{N}
{
  T = new double[d];
  for (int i = 0; i < d; i++)
    T[i] = 2.0*t*(i+1);
}

inline void Dos::setDos()
{
  ep.makeEqDist(N, -T[d-1], T[d-1]);
  double fac = 1./(T[0]*std::pow(pi,d));
  if (d == 1) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
      double K = Maths::min(one, sqr(ep[i]/T[0]));
      D[i] = fac/sqrt(1.0-K);
    }
  } else if (d == 2) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
      double K = Maths::min(one, sqrt((T[1]-ep[i])*(T[1]+ep[i]))/T[1]);
      D[i] = fac*Ellip::complete1(K);
    }
  } else if (d == 3) {
    int j = 0;
    int w0 = ep.fInd(0.0, j);
    int w1 = ep.fInd(T[0],j)+1;
    int w2 = ep.fInd(T[2],j);

    #pragma omp parallel for
    for (int i = w0; i < w1; i++) D[i] = fac*integral0(0.0, pi, ep[i]);
    #pragma omp parallel for
    for (int i = w1; i < w2; i++) D[i] = fac*integral0(0.0, acos((ep[i]-T[1])/T[0]), ep[i]);
    #pragma omp parallel for
    for (int i = 0 ; i < w0; i++) D[i] = D[N-1-i];
  } else {
    cerr << "onle dimension 1, 2, 3 is implemented." << endl;
    exit(1);
  }
}

inline double Dos::RealG0(double om_)
{
  double fac = 1./(T[0]*std::pow(pi,d-1));
  if (d == 1) {
  } else if (d == 2) {
  } else if (d == 3) {
    fac *= sgn(om_);
    double om = fabs(om_);
    
    if (om < T[0]) {
      double ac = acos(om/T[0]);
      return fac*(-integral2(0.0, ac, om) + integral1(ac, pi, om));
    } else if (om < T[2]) {
      double ac = acos((om-T[1])/T[0]);
      return fac*(integral1(0.0, ac, om) + integral3(ac, pi, om));
    } else {
      return fac*(integral3(0.0, pi, om));
    }
  } else {
    cerr << "onle dimension 1, 2, 3 is implemented." << endl;
    exit(1);
  }
  return 0;
}

inline double Dos::ImagG0(double om_)
{
  double om = fabs(om_);
  double fac = -1./(T[0]*std::pow(pi,d-1));
  if (d == 1) {
    double K = Maths::min(one, sqr(om/T[0]));
    return fac/sqrt(1.0-K);
  } else if (d == 2) {
    double K = Maths::min(one, sqrt((T[1]-om)*(T[1]+om))/T[1]);
    return fac*Ellip::complete1(K);
  } else if (d == 3) {
    if (om < T[0]) {
      return fac*integral0(0.0, pi, om);
    } else if (om < T[2]) {
      return fac*integral0(0.0, acos((om-T[1])/T[0]), om);
    } else {
      return 0.0;
    }
  } else {
    cerr << "onle dimension 1, 2, 3 is implemented." << endl;
    exit(1);
  }
}

inline double Dos::integral0(double a, double b, double e)
{ // \int_a^b d\phi K( \frac{sqrt{k^2 - 1}}{k} ), k = \frac{4t}{\e - 2t\cos\phi}
  // Note: K(1) happens when \phi_c = cos^{-1}(\e / 2t).
  //       This occurs if \e is in [0,2t].
  M1d phi;
  phi.makeEqDist(n, a, b);
  double sum = 0.0;
  for (int j = 0; j < phi.size(); j++) {
    double k = T[1]/(e - T[0]*cos(phi[j]));
    double K = Maths::min(one, fabs(sqrt(k*k-1)/k));
    sum += phi.dx(j)*Ellip::complete1(K);
  }
  return sum;
}

inline double Dos::integral1(double a, double b, double e)
{ // \int_a^b d\phi K( \frac{1}{k} ), k = \frac{4t}{\e - 2t\cos\phi}
  M1d phi;
  phi.makeEqDist(n, a, b);
  double sum = 0.0;
  for (int j = 0; j < phi.size(); j++) {
    double k = T[1]/(e - T[0]*cos(phi[j]));
    double K = Maths::min(one, 1.0/k);
    sum += phi.dx(j)*Ellip::complete1(K);
  }
  return sum;
}
inline double Dos::integral2(double a, double b, double e)
{ // \int_a^b d\phi K( \frac{1}{|k|} ), k = \frac{4t}{\e - 2t\cos\phi}
  M1d phi;
  phi.makeEqDist(n, a, b);
  double sum = 0.0;
  for (int j = 0; j < phi.size(); j++) {
    double k = T[1]/(e - T[0]*cos(phi[j]));
    double K = Maths::min(one, 1.0/fabs(k));
    sum += phi.dx(j)*Ellip::complete1(K);
  }
  return sum;
}

inline double Dos::integral3(double a, double b, double e)
{ // \int_a^b d\phi k K( k ), k = \frac{4t}{\e - 2t\cos\phi}
  M1d phi;
  phi.makeEqDist(n, a, b);
  double sum = 0.0;
  for (int j = 0; j < phi.size(); j++) {
    double k = T[1]/(e - T[0]*cos(phi[j]));
    double K = Maths::min(one, k);
    sum += phi.dx(j)*k*Ellip::complete1(K);
  }
  return sum;
}

#endif /* end of include guard: NONINTDOS_H */
