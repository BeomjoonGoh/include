#ifndef NONINTDOS_H
#define NONINTDOS_H

#include <iostream>
#include <algorithm>

#include "maths.h"
#include "function.h"
#include "mesh.h"
#include "specialfunctions/elliptic.h"
#include "derule.h"

#include <omp.h>

using namespace Maths;

class Dos
{ // Non interacting density of states of simple cupic lattice for 1, 2, and 3 dimension.
  // E. N. Economou, Green's Functions in Quantum Physics, 3rd Ed., (Springer, 2006), p. 87.
  // T. Morita and T. Horiguchi, J. Math. Phys. 12 981 (1971).
  private:
    int d;
    int N;
    double *T;
    F1d<double> D;
    int w0, w1;
    const double one;

  public:
    M1d ep;

  public:
    Dos() : T{nullptr}, one{1.0-absEpsilon} { }
    Dos(int d_, int N_, double t) : T{nullptr}, one{1.0-absEpsilon} { init(d_, N_, t); }
    Dos(const Dos &dos);
    ~Dos() { delete[] T; T = nullptr;}

    Dos& operator= (const Dos &dos);
    double& operator[] (const int i) { return D[i]; }
    const double& operator[] (const int i) const { return D[i]; }
    double de(int i) const { return ep.dx(i); }
    double d1e(int i) const { return ep.d1x(i); }

    void init(int d_, int N_, double t);
    int size() const { return N; }
    void setDos();
    double RealG0(double om);
    double ImagG0(double om);
    void print(const std::string &outputf, const std::string &comment = "") { D.print(ep,outputf,comment); } 

};

// {{{
class Func_Base
{ // Base class for functor classes Func0..3. In all cases, k = \frac{4t}{\e - 2t\cos\phi}, l = 1/k.
  // The functor classes are for quadrature d = 3 case.
  public:
    const double eps = 1e-9;
  protected:
    const double *T;
    const double one;
    const double e;

    Func_Base(double *T_, double one_, double e_) : T{T_}, one{one_}, e{e_} { }
    double k(double e, double phi) { return T[1]/(e - T[0]*cos(phi)); }
    double l(double e, double phi) { return (e - T[0]*cos(phi))/T[1]; }
};

class Func0 : public Func_Base
{ // K( \frac{\sqrt{k^2 - 1}}{k} ) = K(\sqrt{1-l^2})
  // Note: K(1) happens if \e \in [0,2t] at \phi_c = \cos^{-1}(\e / 2t)
  public:
    Func0(double *T_, double one_, double e_) : Func_Base{T_, one_, e_} { }
    double operator()(double phi, double dum) { return Ellip::complete1( Maths::min(one, sqrt(1-sqr(l(e, phi))))); }
};

class Func1 : public Func_Base
{ // K( \frac{1}{k} ) = K(l)
  // Note: K(1) happens if \e \in [2t,6t] at \phi_c = \cos^{-1}( \frac{\e - 4t}{2t} )
  public:
    Func1(double *T_, double one_, double e_) : Func_Base{T_, one_, e_} { }
    double operator()(double phi, double dum) { return Ellip::complete1( Maths::min(one, l(e,phi)) ); }
};

class Func2 : public Func_Base
{ // K( \frac{1}{|k|} ) = K(|l|)
  // Note: K(1) happens if \e \in [ 2t, 6t] at \phi_c = \cos^{-1}( \frac{\e - 4t}{2t} )
  //                    if \e \in [-6t,-2t] at \phi_c = \cos^{-1}( \frac{\e + 4t}{2t} )
  public:
    Func2(double *T_, double one_, double e_) : Func_Base{T_, one_, e_} { }
    double operator()(double phi, double dum) { return Ellip::complete1( Maths::min(one, abs(l(e,phi))) ); }
};

class Func3 : public Func_Base
{ // k K( k )
  // Note: K(1) happens if \e \in [2t,6t] at \phi_c = \cos^{-1}( \frac{\e - 4t}{2t} )
  public:
    Func3(double *T_, double one_, double e_) : Func_Base{T_, one_, e_} { }
    double operator()(double phi, double dum) { double K = k(e,phi); return K*Ellip::complete1( Maths::min(one, K) ); }
};
//}}}

Dos::Dos(const Dos &dos) : d(dos.d), N(dos.N), T(dos.d ? new double[dos.d] : nullptr), D(dos.D),
                           w0(dos.w0), w1(dos.w1), one(dos.one), ep(dos.ep)
{
  if (T) 
    std::copy(dos.T, dos.T+dos.d, T);
}

Dos& Dos::operator= (const Dos &dos)
{
  if (this == &dos)
    return *this;
  this->d = dos.d;
  this->N = dos.N;

  if (this->T) delete[] this->T;
  this->T = new double[dos.d];
  std::copy(dos.T, dos.T+dos.d, this->T);

  this->D = dos.D;
  this->w0 = dos.w0;
  this->w1 = dos.w1;
  this->ep = dos.ep;

  return *this;
}

void Dos::init(int d_, int N_, double t)
{
  d = d_;
  N = N_>>1<<1;
  D.resize(N);
  T = new double[d];
  for (int i = 0; i < d; i++)
    T[i] = 2.0*t*(i+1); // 2t, 4t, 6t, ... 
}

inline void Dos::setDos()
{
  double fac = 1./(T[0]*std::pow(pi,d));
  switch (d) {
    case 1: {
      double dE = 2*T[d-1]/(N-1);
      ep.makeEqDist(N, -T[d-1]+dE, T[d-1]-dE);
      #pragma omp parallel for
      for (int i = 0; i < N; i++) D[i] = fac/sqrt(1.0-Maths::min(one, sqr(ep[i]/T[0])));
      break;
    }
    case 2: {
      ep.makeEqDist(N, -T[d-1], T[d-1]);
      #pragma omp parallel for
      for (int i = 0; i < N; i++) D[i] = fac*Ellip::complete1(Maths::min(one, sqrt((T[1]-ep[i])*(T[1]+ep[i]))/T[1]));
      break;
    }
    case 3: {
      ep.makeEqDist(N, -T[d-1], T[d-1]);

      int j = 0;
      w0 = ep.fInd(0.0, j);
      w1 = ep.fInd(T[0],j)+1;

      #pragma omp parallel for
      for (int i = w0; i < w1; i++) {
        Func0 f0(T,one,ep[i]);
        D[i] = fac*(integral(0.0, acos(ep[i]/T[0]), f0, f0.eps) + integral(acos(ep[i]/T[0]), pi, f0, f0.eps));
      }
      #pragma omp parallel for
      for (int i = w1; i < N; i++) {
        Func0 f0(T,one,ep[i]);
        D[i] = fac*integral(0.0, acos((ep[i]-T[1])/T[0]), f0, f0.eps);
      }
      #pragma omp parallel for
      for (int i = 0 ; i < w0; i++) D[i] = D[N-1-i];
      break;
    }
    default: {
      quitif(true, "only for dimension 1, 2, 3 are implemented.");
    }
  }
}

inline double Dos::RealG0(double om_)
{
  double fac = sgn(om_)/(T[0]*std::pow(pi,d-1));
  double om = abs(om_);
  switch (d) {
    case 1: {
      if (om < T[0])
        return 0.0;
      return fac/sqrt(sqr(om/T[0])-one);
    }
    case 2: {
      if (om < T[0]) {
        double l = om/T[0];
        return fac*Ellip::complete1(Maths::min(one,l));
      } else {
        double k = T[0]/om;
        return fac*k*Ellip::complete1(Maths::min(one,k));
      }
    }
    case 3: {
      Func1 f1(T,one,om);
      Func2 f2(T,one,om);
      Func3 f3(T,one,om);
      if (om < T[0]) {
        double ac = acos(om/T[0]);
        return fac*(-integral(0.0, ac, f2, f2.eps) + integral(ac, pi, f1, f1.eps));
      } else if (om < T[2]) {
        double ac = acos((om-T[1])/T[0]);
        return fac*(integral(0.0, ac, f1, f1.eps) + integral(ac, pi, f3, f3.eps));
      } else {
        return fac*(integral(0.0, pi, f3, f3.eps));
      }
    }
    default: {
      quitif(true, "only for dimension 1, 2, 3 are implemented.");
    }
  }
}

inline double Dos::ImagG0(double om_)
{
  double om = abs(om_);
  double fac = -1./(T[0]*std::pow(pi,d-1));
  switch (d) {
    case 1: {
      if (om < T[0])
        return fac/sqrt(1.0-Maths::min(one, sqr(om/T[0])));
      return 0.0;
    }
    case 2: {
      if (om < T[0])
        return fac*Ellip::complete1(Maths::min(one, sqrt((T[1]-om)*(T[1]+om))/T[1]));
      return 0.0;
    }
    case 3: {
      int j0 = w1;
      if (om <= ep[w1-1])
        return -pi*D.line(ep.getInterp(om, j0, &M1d::finD));
      j0 = w1-2;
      if (om <= T[0])
        return -pi*D.line(Interp{j0, (om-ep[j0])*ep.d1x(j0)});
      if (om < ep[w1]) { 
        double DT0 = D.line(Interp{j0, (T[0]-ep[j0])*ep.d1x(j0)});
        return -pi*(DT0 - (D[w1]-DT0)*(om-T[0])/(ep[w1]-T[0]));
      }
      if (om <= T[2])
        return -pi*D.line(ep.getInterp(om, j0, &M1d::fInd));
      return 0.0;
    }
    default: {
      quitif(true, "only for dimension 1, 2, 3 are implemented.");
    }
  }
}

/*
void DMFT::readDos(string inputf)
{
  F2d<double> Dos;
  M1d ep;
  Dos.resize(Nb);
  if (!Dos.readx(ep,inputf,"non-interacting D density of states")) {
    TightBinding tb{Param::dim, Param::latUnitVec, Param::kmesh};
    tb.computeHk(Nb, Param::hopT, Param::orbitals);
    kSize = tb.size();
    weight = 1./tb.size();

    ek.resize(Nb, Nb);
    for (int b = 0; b < Nb; b++) {
      for (int c = 0; c < Nb; c++) {
        ek(b,c).resize(kSize);
        for (int k = 0; k < kSize; k++)
          ek(b,c)[k] = tb.Hk[b][c][k].real();
      }
    }
    clog << "Making density of states" << endl;
    F1d<double> Bmax{Nb}; F1d<double> Bmin{Nb};
    for (int b = 0; b < Nb; b++) {
      Bmax[b] = *max_element(ek(b,b).begin(),ek(b,b).end());
      Bmin[b] = *min_element(ek(b,b).begin(),ek(b,b).end());
    }
    int m0 = om.fInd(*min_element(Bmin.begin(),Bmin.end()))-10;
    int m1 = om.fInd(*max_element(Bmax.begin(),Bmax.end()))+11;
    Dos.resize(Nb, m1-m0+1);
    for (int b = 0; b < Nb; b++) {
      double h2 = sqr(Param::eta);
      #pragma omp parallel for
      for (int i = m0; i <= m1; i++) {
        double sum = 0;
        for (int k = 0; k < kSize; k++) {
          sum += Param::eta/(sqr(om[i]-ek(b,b)[k]) + h2);
        }
        Dos(b,i-m0) = sum*weight/pi;
      }
    }
    ep.resize(Dos.size());
    for (int e = 0; e < ep.size(); e++)
      ep[e] = om[e+m0];
    ep.Set();

    Dos.print(ep,inputf,Param::HamInfo());
  }
}
*/
#endif /* end of include guard: NONINTDOS_H */
