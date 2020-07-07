#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include "assert.h"
#include "function.h"
#include "mesh.h"

class Coeff
{
  public:
    int j;
    double Cn, Co;
    double xp1, xm1;

    int i; 
    double d, p1, p2, q1, q2;
    bool linear;
  
    Coeff() { }
    ~Coeff() { }

    void init(const double z, const M1d &x, int &ii);
    void calculate_0(const double z, const M1d &x, int &ii);
    void calculate_j(const double z, const int j_, const M1d &x, int &ii);
    void calculate_N(const double z, const M1d &x);

  private:
    void Setpq(const double xp1, const int i, const M1d &x);
};

template <typename T>
class Moments
{
  private:
    F1d<T> F0;
    F1d<T> F1;
    F1d<T> F2;

    static const int N = 2;
    T Fxz[N];
    T dF_old[N];
    T dF_new[N];

  public:
    Moments() { }
    Moments(const F1d<T> &f, const M1d &x);
    ~Moments() { }

    void Set(const F1d<T> &f, const M1d &x);
    void init(const Coeff &c);
    T getAxdx_j(const Coeff &c, const M1d &x);
    T getAxdx_N(const Coeff &c, const M1d &x);

  private:
    void resize(int N = 0);
    void quadInterp(const Coeff &c, T *F);
    void update(const T *F);
};

namespace Conv
{
  template <typename T> void makeAverage(double z, const M1d &om, Moments<T> &M, F1d<T> &A);
  template <typename T> T convolute(const F1d<T> &g, const F1d<T> &A);
  template <typename T> T convolute(const F1d<T> &A, const F1d<T> &B, const M1d &x, const int i, const int m);
}

// Coeff
void Coeff::init(const double z, const M1d &x, int &ii)
{
  xp1 = x[0]-z;
  Interp I = x.getInterp(xp1, ii, &M1d::fInd);
  i = I.j;
  d = I.d;

  Setpq(xp1,i,x);
}

inline void Coeff::calculate_0(const double z, const M1d &x, int &ii)
{
  j = 0;
  Cn = x.d1x(j);
  Co = 0;
  xp1 = x[j+1]-z;
  xm1 = 0;

  Interp I = x.getInterp(xp1, ii, &M1d::fInd);
  i = I.j;
  d = I.d;

  Setpq(xp1,i,x);
  linear = false;
}

inline void Coeff::calculate_j(const double z, const int j_, const M1d &x, int &ii)
{
  const double small = 0.01;
  j = j_;
  Cn = x.d1x(j);
  Co = x.d1x(j-1);
  xp1 = x[j+1]-z;
  xm1 = x[j-1]-z;

  Interp I = x.getInterp(xp1, ii, &M1d::fInd);
  i = I.j;
  d = I.d;

  Setpq(xp1,i,x);
  linear = x.dx(j) < small*x.dx(i);
}

inline void Coeff::calculate_N(const double z, const M1d &x)
{
  j = x.size()-1;
  Cn = 0;
  Co = x.d1x(j-1);
  xp1 = 0;
  xm1 = x[j-1]-z;
  linear = false;
}

inline void Coeff::Setpq(const double xp1, const int i, const M1d &x)
{
  double dx = xp1 - x[i];
  double dx2= dx*dx;
  p1 = x.d1x(i)*dx2/2.0;
  q1 = dx - p1;

  p2 = x.d1x(i)*dx2*(2.0*xp1 + x[i])/6.0;
  q2 = dx*(xp1+x[i])/2.0 - p2;
}

// Moments
template <typename T>
inline Moments<T>::Moments(const F1d<T> &f, const M1d &x)
{
  Set(f,x);
}

template <typename T>
inline void Moments<T>::Set(const F1d<T> &f, const M1d& x)
{
  resize(x.size());
  assert0(f.size() == x.size());
  F0[0] = f[0];
  F1[0] = 0;
  F2[0] = 0;
  for (int i = 1; i < x.size(); i++) {
    double dx = 1./x.d1x(i-1);
    F0[i] = f[i];
    F1[i] = F1[i-1] + (F0[i] + F0[i-1])*dx/2.0;
    F2[i] = F2[i-1] + (F0[i]*(2.0*x[i]+x[i-1]) + F0[i-1]*(x[i]+2.0*x[i-1]))*dx/6.0;
  }
}

template <typename T>
inline void Moments<T>::init(const Coeff &c)
{
  quadInterp(c, Fxz);
  for (int k = 0; k < N; k++)
    dF_new[k] = Fxz[k];
}

template <typename T>
inline T Moments<T>::getAxdx_j(const Coeff &c, const M1d &x)
{
  T F[N];
  quadInterp(c, F);
  update(F);

  if (c.linear) return F0.line(Interp(c.i,c.d)) * x.dx(c.j);
  else          return c.Co*(dF_old[1] - c.xm1*dF_old[0]) - c.Cn*(dF_new[1] - c.xp1*dF_new[0]);
}

template <typename T>
inline T Moments<T>::getAxdx_N(const Coeff &c, const M1d &x)
{
  return c.Co*(dF_old[1] - c.xm1*dF_old[0]);
}
  
template <typename T>
inline void Moments<T>::resize(int N)
{
  F0.resize(N);
  F1.resize(N);
  F2.resize(N);
}

template <typename T>
inline void Moments<T>::quadInterp(const Coeff &c, T *F)
{
  F[0] = F1[c.i] + c.q1*F0[c.i] + c.p1*F0[c.i+1];
  F[1] = F2[c.i] + c.q2*F0[c.i] + c.p2*F0[c.i+1];
}

template <typename T>
inline void Moments<T>::update(const T *F)
{
  for (int k = 0; k < N; k++) {
    dF_old[k] = dF_new[k];
    dF_new[k] = F[k] - Fxz[k];
    Fxz[k] = F[k];
  }
}

// Conv
template <typename T>
void Conv::makeAverage(double z, const M1d &om, Moments<T> &M, F1d<T> &A)
{
  int ii = 0;
  Coeff c;
  c.init(z,om,ii);
  M.init(c);
  c.calculate_0(z, om, ii);
  A[0] = M.getAxdx_j(c, om);
  for (int j = 1; j < om.size()-1; j++) {
    c.calculate_j(z, j, om, ii);
    A[j] = M.getAxdx_j(c, om);
  }
  c.calculate_N(z, om);
  A[om.size()-1] = M.getAxdx_N(c, om);
}

template <typename T>
inline T Conv::convolute(const F1d<T> &g, const F1d<T> &A)
{ // Instead, use BLAS DGEMM with F2d g, F2d A?
  // or use BLAS DDOT with F1d g, F1d A? 
  // or std::inner_product(F1d g, F1d A, 0.0)?
  assert0(g.size() == A.size());
  T sum(0);
  for (int i = 0; i < g.size(); i++)
    sum += g[i]*A[i];
  return sum;
}

template <typename T>
inline T Conv::convolute(const F1d<T> &A, const F1d<T> &B, const M1d &x, const int i, const int m)
{ // returns value of: C(w_i) = int dx A(x)B(x-w_i) for equidistance mesh
  T sum(0);
  for (int j = 0; j < x.size(); j++) {
    int ind = j-i+(m-1);
    if (ind < 0)          continue;
    if (ind > x.size()-1) break;
    sum += x.dx(j)*A[j]*B[ind];
  }
  return sum;
}

#endif /* end of include guard: CONVOLUTION_H */
