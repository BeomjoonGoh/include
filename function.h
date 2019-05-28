#ifndef FUNCTION_H
#define FUNCTION_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <functional>

#include "util.h"
#include "maths.h"
#include "complex.h"
#include "assert.h"

template <typename T>
class F1d
{
  protected:
    int N;              // number of elements
    T *v;               // function
    // Since the class is templated, the usage of member variables in the derived class must be like: this->N

  public:
    F1d() : N{0}, v{nullptr} { };
    F1d(int N_);
    F1d(int N_, T vi);
    F1d(const F1d<T> &f);
    virtual ~F1d() { delete[] v; v = nullptr;}

    // operators {{{
    T& operator[] (const int i) { assert(0<=i && i<N, "i="<<i<<" N="<< N); return v[i]; }
    const T& operator[] (const int i) const { assert(0<=i && i<N, "i="<<i<<" N="<<N); return v[i]; }

    F1d<T> operator- () const;
    F1d<T>& operator=  (const F1d<T> &f);
    F1d<T>& operator+= (const F1d<T> &f);
    F1d<T>& operator+= (const T &a);
    F1d<T>& operator-= (const F1d<T> &f);
    F1d<T>& operator-= (const T &a);
    F1d<T>& operator*= (const F1d<T> &f);
    F1d<T>& operator*= (const T &a);
    F1d<T>& operator/= (const F1d<T> &f);
    F1d<T>& operator/= (const T &a);

    template <typename U> friend F1d<U> operator+ (const F1d<U> &f1, const F1d<U> &f2);
    template <typename U> friend F1d<U> operator+ (const F1d<U> &f, const U &a);
    template <typename U> friend F1d<U> operator+ (const U &a, const F1d<U> &f);
    template <typename U> friend F1d<U> operator- (const F1d<U> &f1, const F1d<U> &f2);
    template <typename U> friend F1d<U> operator- (const F1d<U> &f, const U &a);
    template <typename U> friend F1d<U> operator- (const U &a, const F1d<U> &f);
    template <typename U> friend F1d<U> operator* (const F1d<U> &f1, const F1d<U> &f2);
    template <typename U> friend F1d<U> operator* (const F1d<U> &f, const U &a);
    template <typename U> friend F1d<U> operator* (const U &a, const F1d<U> &f);
    template <typename U> friend F1d<U> operator/ (const F1d<U> &f1, const F1d<U> &f2);
    template <typename U> friend F1d<U> operator/ (const F1d<U> &f, const U &a);
    template <typename U> friend F1d<U> operator/ (const U &a, const F1d<U> &f);
    // operators }}}

    T* begin() const {return v;}
    T* end() const {return v+N;}

    T& First() { return v[0]; }
    T& Last() { return v[N-1]; }
    const T& First() const { return v[0]; }
    const T& Last() const { return v[N-1]; }

    int size() const { return N; }
    void init(int N_, T* v_) { this->N = N_; this->v = v_; }
    void resize(int N_ = 0);
    void Set(int N_, T vi);
    void push_back(const T &last);
    T Sum();
    double Rms();
    template <typename U> friend double normalizedDiff(const F1d<U> &f1, const F1d<U> &f2);

    template <typename M1d> void Print(const M1d &x, const std::string &outputf, const std::string &comment = "");
    void Print(const std::string &outputf, const std::string &comment = "");
    bool Read(const std::string &inputf, const std::string &descrip = "");
    template <typename M1d> bool Readx(M1d &x, const std::string &inputf, const std::string &descrip = "");

    template <typename U, typename M1d>
    friend F1d<U> kramerskronig(const M1d &x, const F1d<U> &f);
    template <typename U> friend U integralTZ (F1d<U> &f, double dx, int i, int j);
    template <typename U> friend U integralTZ2(F1d<U> &f, double dx, int i, int j);
    template <typename M1d> void cspline(const M1d &X1, const F1d &Y1, const M1d &X2);
    template <typename Interp> T line(const Interp &jd) const;

  private:
    template <typename M1d> void getYD_gen(const M1d &X, const F1d<T> &Y, F1d<T> &YD);
    void tridiag_gen(const F1d<T> &A, const F1d<T> &B, const F1d<T> &C, F1d<T> &D);
};

template <typename T>
class F2d
{
  protected:
    void *mem;
    int Nd;             // number of elements
    int N;              // number of elements
    F1d<T> *v;          // functions

  public:
    F2d() : mem{nullptr}, Nd{0}, N{0}, v{nullptr} { }
    F2d(int Nd_, int N_);
    F2d(int Nd_, int N_, T vi);
    F2d(int Nd_, const F1d<T> &f);
    F2d(const F2d<T> &F);
    virtual ~F2d() { operator delete(mem); mem = nullptr; }

    // operators {{{
    F1d<T>& operator[] (const int i) { assert(0<=i && i<Nd, "i="<<i<< " Nd="<<Nd); return v[i]; }
    const F1d<T>& operator[] (const int i) const { assert(0<=i && i<Nd, "i="<<i<< " Nd="<<Nd); return v[i]; }
    T& operator()(const int i, const int j) { assert(0<=i && i<Nd && 0<=j && j<N,"i,j="<<i<<","<<j<<" Nd,N="<<Nd<<","<<N); return v[i][j]; }
    const T& operator()(const int i, const int j) const { assert(0<=i && i<Nd && 0<=j && j<N,"i,j="<<i<<","<<j<<" Nd,N="<<Nd<<","<<N); return v[i][j]; }

    F2d<T> operator- () const;
    F2d<T>& operator=  (const F2d<T> &F);
    F2d<T>& operator+= (const F2d<T> &F);
    F2d<T>& operator+= (const F1d<T> &f);
    F2d<T>& operator+= (const T &a);
    F2d<T>& operator-= (const F2d<T> &F);
    F2d<T>& operator-= (const F1d<T> &f);
    F2d<T>& operator-= (const T &a);
    F2d<T>& operator*= (const F2d<T> &F);
    F2d<T>& operator*= (const F1d<T> &f);
    F2d<T>& operator*= (const T &a);
    F2d<T>& operator/= (const F2d<T> &F);
    F2d<T>& operator/= (const F1d<T> &f);
    F2d<T>& operator/= (const T &a);

    template <typename U> friend F2d<U> operator+ (const F2d<U> &F1, const F2d<U> &F2);
    template <typename U> friend F2d<U> operator+ (const F2d<U> &F1, const F1d<U> &f2);
    template <typename U> friend F2d<U> operator+ (const F1d<U> &f1, const F2d<U> &F2);
    template <typename U> friend F2d<U> operator+ (const F2d<U> &F, const U &a);
    template <typename U> friend F2d<U> operator+ (const U &a, const F2d<U> &F);
    template <typename U> friend F2d<U> operator- (const F2d<U> &F1, const F2d<U> &F2);
    template <typename U> friend F2d<U> operator- (const F2d<U> &F1, const F1d<U> &f2);
    template <typename U> friend F2d<U> operator- (const F1d<U> &f1, const F2d<U> &F2);
    template <typename U> friend F2d<U> operator- (const F2d<U> &F, const U &a);
    template <typename U> friend F2d<U> operator- (const U &a, const F2d<U> &F);
    template <typename U> friend F2d<U> operator* (const F2d<U> &F1, const F2d<U> &F2);
    template <typename U> friend F2d<U> operator* (const F2d<U> &F1, const F1d<U> &f2);
    template <typename U> friend F2d<U> operator* (const F1d<U> &f1, const F2d<U> &F2);
    template <typename U> friend F2d<U> operator* (const F2d<U> &F, const U &a);
    template <typename U> friend F2d<U> operator* (const U &a, const F2d<U> &F);
    template <typename U> friend F2d<U> operator/ (const F2d<U> &F1, const F2d<U> &F2);
    template <typename U> friend F2d<U> operator/ (const F2d<U> &F1, const F1d<U> &f2);
    template <typename U> friend F2d<U> operator/ (const F1d<U> &f1, const F2d<U> &F2);
    template <typename U> friend F2d<U> operator/ (const F2d<U> &F, const U &a);
    template <typename U> friend F2d<U> operator/ (const U &a, const F2d<U> &F);
    // operators }}}

    T* begin() {return (*v).begin();}

    int size() const { return N; }
    int sizeD() const { return Nd; }
    void resize(int Nd_ = 0, int N_ = 0); 
    void Set(int Nd_, int N_, T vi);
    double Rms();
    template <typename U> friend double normalizedDiff(const F2d<U> &F1, const F2d<U> &F2);

    template <typename M1d> void Print(const M1d &x, const std::string &outputf, const std::string &comment = "");
    bool Read(const std::string &inputf, const std::string &descrip = "");
    bool Readr(const std::string &inputf);
    template <typename M1d> bool Readx(M1d &x, const std::string &inputf, const std::string &descrip = "");

    template <typename U, typename M1d> friend F2d<U>
    kramerskronig(const M1d &x, const F2d<U> &F);
};

// ============  F1d  ============
// ======  Constructors
template <typename T>
inline F1d<T>::F1d(int N_)
{
  assert0(N_ >= 0);
  N = N_;
  v = new T[N];
}
template <typename T>
inline F1d<T>::F1d(int N_, T vi)
{
  assert0(N_ >= 0);
  N = N_;
  v = new T[N];
  for (int i = 0; i < N; i++)
    v[i] = vi;
}
template <typename T>
inline F1d<T>::F1d(const F1d<T> &f)
{
  N = f.N;
  v = new T[N];
  std::copy(f.v, f.v+f.N, v);
}

// ====== Operator overloadings
template <typename T>
inline F1d<T> F1d<T>::operator- () const
{
  F1d<T> f{this->N};
  for (int i = 0; i < f.N; i++)
    f.v[i] = -static_cast<T>(1)*this->v[i];
  return f;
}
template <typename T>
inline F1d<T>& F1d<T>::operator= (const F1d<T> &f)
{
  if (this == &f)
    return *this;
  resize(f.N);
  std::copy(f.v, f.v+f.N, this->v);
  return *this;
}
template <typename T>
inline F1d<T>& F1d<T>::operator+= (const F1d<T> &f)
{
  assert0(this->N == f.N);
  for (int i = 0; i < this->N; i++)
    this->v[i] += f.v[i];
  return *this;
} 
template <typename T>
inline F1d<T>& F1d<T>::operator+= (const T &a)
{
  for (int i = 0; i < this->N; i++)
    this->v[i] += a;
  return *this;
} 
template <typename T>
inline F1d<T>& F1d<T>::operator-= (const F1d<T> &f)
{
  assert0(this->N == f.N);
  for (int i = 0; i < this->N; i++)
    this->v[i] -= f.v[i];
  return *this;
} 
template <typename T>
inline F1d<T>& F1d<T>::operator-= (const T &a)
{
  for (int i = 0; i < this->N; i++)
    this->v[i] -= a;
  return *this;
} 
template <typename T>
inline F1d<T>& F1d<T>::operator*= (const F1d<T> &f)
{
  assert0(this->N == f.N);
  for (int i = 0; i < this->N; i++)
    this->v[i] *= f.v[i];
  return *this;
} 
template <typename T>
inline F1d<T>& F1d<T>::operator*= (const T &a)
{
  for (int i = 0; i < this->N; i++)
    this->v[i] *= a;
  return *this;
} 
template <typename T>
inline F1d<T>& F1d<T>::operator/= (const F1d<T> &f)
{
  assert0(this->N == f.N);
  for (int i = 0; i < this->N; i++) {
    assert(f.v[i] != 0., "Divide by zero");
    this->v[i] /= f.v[i];
  }
  return *this;
} 
template <typename T>
inline F1d<T>& F1d<T>::operator/= (const T &a)
{
  assert(a != 0., "Divide by zero");
  for (int i = 0; i < this->N; i++)
    this->v[i] /= a;
  return *this;
} 

template <typename T>
inline F1d<T> operator+ (const F1d<T> &f1, const F1d<T> &f2)
{
  assert0(f1.N == f2.N);
  F1d<T> g{f1.N};
  for (int i = 0; i < f1.N; i++)
    g[i] = f1[i] + f2[i];
  return g;
} 
template <typename T>
inline F1d<T> operator+ (const F1d<T> &f, const T &a)
{
  F1d<T> g{f.N};
  for (int i = 0; i < f.N; i++)
    g[i] = f[i] + a;
  return g;
} 
template <typename T>
inline F1d<T> operator+ (const T &a, const F1d<T> &f)
{
  return f + a;
} 

template <typename T>
inline F1d<T> operator- (const F1d<T> &f1, const F1d<T> &f2)
{
  assert0(f1.N == f2.N);
  F1d<T> g{f1.N};
  for (int i = 0; i < f1.N; i++)
    g[i] = f1[i] - f2[i];
  return g;
} 
template <typename T>
inline F1d<T> operator- (const F1d<T> &f, const T &a)
{
  return f + (-a);
}
template <typename T>
inline F1d<T> operator- (const T &a, const F1d<T> &f) {
  F1d<T> g{f.N};
  for (int i = 0; i < f.N; i++)
    g[i] = a - f[i];
  return g;
} 

template <typename T>
inline F1d<T> operator* (const F1d<T> &f1, const F1d<T> &f2)
{
  assert0(f1.N == f2.N);
  F1d<T> g{f1.N};
  for (int i = 0; i < f1.N; i++)
    g[i] = f1[i] * f2[i];
  return g;
} 
template <typename T>
inline F1d<T> operator* (const F1d<T> &f, const T &a)
{
  F1d<T> g{f.N};
  for (int i = 0; i < f.N; i++)
    g[i] = f[i] * a;
  return g;
} 
template <typename T>
inline F1d<T> operator* (const T &a, const F1d<T> &f)
{
  return f * a;
} 

template <typename T>
inline F1d<T> operator/ (const F1d<T> &f1, const F1d<T> &f2)
{
  assert0(f1.N == f2.N);
  F1d<T> g{f1.N};
  for (int i = 0; i < f1.N; i++) {
    assert(f2[i] != 0., "Divide by zero");
    g[i] = f1[i] / f2[i];
  }
  return g;
}
template <typename T>
inline F1d<T> operator/ (const F1d<T> &f, const T &a)
{
  assert(a != 0., "Divide by zero");
  return f * (1./a);
} 
template <typename T>
inline F1d<T> operator/ (const T &a, const F1d<T> &f)
{
  F1d<T> g{f.N};
  for (int i = 0; i < f.N; i++) {
    assert(f[i] != 0., "Divide by zero");
    g[i] = a / f[i];
  }
  return g;
}

template <typename T>
inline void F1d<T>::resize(int N_)
{
  assert0(N_ >= 0);
  if (N_ > N) {
    if (v) delete[] v;
    v = new T[N_];
  }
  N = N_;
}
template <typename T>
inline void F1d<T>::Set(int N_, T vi)
{
  resize(N_);
  for (int i = 0; i < N; i++)
    v[i] = vi;
}
template <typename T>
inline void F1d<T>::push_back(const T &last)
{
  T *vtmp = new T[N+1];
  std::copy(v,v+N,vtmp);
  vtmp[N] = last;
  delete[] v;

  N++;
  v = new T[N];
  std::copy(vtmp,vtmp+N,v);
  delete[] vtmp;
}
template <typename T>
inline T F1d<T>::Sum()
{
  T sum{0};
  for (int i = 0; i < N; i++)
    sum += v[i];
  return sum;
}
template <typename T>
inline double F1d<T>::Rms()
{
  double sum{0};
  for (int i = 0; i < N; i++)
    sum += (conj(v[i])*v[i]).real();
  return abs(sqrt(sum))/N;
}
template <typename T>
double normalizedDiff(const F1d<T> &f1, const F1d<T> &f2)
{
  assert0(f1.N == f2.N);
  double dif{0.}, nor{0.};
  #pragma omp parallel for default(shared) reduction(+:dif,nor)
  for (int i = 0; i < f1.N; i++) {
    dif += Maths::abs(f1[i]-f2[i]);
    nor += Maths::abs(f2[i]);
  }
  double ndif = dif/nor;
  return ndif;
}
template <>
double normalizedDiff(const F1d<compdb> &f1, const F1d<compdb> &f2)
{
  assert0(f1.N == f2.N);
  double dif{0.}, nor{0.};
  #pragma omp parallel for default(shared) reduction(+:dif,nor)
  for (int i = 0; i < f1.N; i++) {
    dif += abs(f1[i]-f2[i]);
    nor += abs(f2[i]);
  }
  double ndif = dif/nor;
  return ndif;
}
template <typename T>
template <typename M1d>
void F1d<T>::Print(const M1d &x, const std::string &outputf, const std::string &comment) {
  assert0(x.size() == N);
  std::clog << "  => Function printed: " << outputf << "\n";

  std::ofstream outf{outputf};
  if (comment != "")
    outf << "# " << comment << '\n';
  outf << std::fixed << std::setprecision(25);
  for (int i = 0; i < N; i++)
    outf << x[i] << '\t' << v[i] << '\n';
}
template <typename T>
void F1d<T>::Print(const std::string &outputf, const std::string &comment)
{
  std::ofstream outf{outputf};
  if (comment != "")
    outf << "# " << comment << '\n';
  outf << std::fixed << std::setprecision(12);
  for (int i = 0; i < N; i++)
    outf << v[i] << '\n';
}
template <typename T>
bool F1d<T>::Read(const std::string &inputf, const std::string &descrip)
{
  std::clog << "Reading " << descrip << " from " << inputf << ", with " << N << " entries.\n";

  std::ifstream inf{inputf};
  if (!inf.good()) {
    std::clog << "Missing input file! " << "No " << inputf << "\n";
    return false;
  }

  double x;
  int i = -1;
  while (inf && ++i < N) {
    if (inf.peek() == '#') inf.ignore(2000,'\n');
    inf >> x >> v[i];
    inf.ignore(2000,'\n');
  }
  if (i!=N) {
    std::clog << "Reading failed: " << inputf << " i(=" << i << ") != N (" << N << ")\n";
    return false;
  }
  Util::getComment(inf,inputf);

  return true;
}
template <typename T>
template <typename M1d> 
bool F1d<T>::Readx(M1d &x, const std::string &inputf, const std::string &descrip)
{
  std::clog << "Reading " << descrip << "with mesh from " << inputf << "\n";

  std::ifstream inf{inputf};
  if (!inf.good()) {
    std::clog << "Missing input file! " << "No " << inputf << "\n";
    return false;
  }

  this->resize(0);
  double xi;
  T vi;
  while (inf) {
    if (inf.peek() == '#') inf.ignore(2000,'\n');
    inf >> xi >> vi;
    this->push_back(vi);
    inf.ignore(2000,'\n');
  }
  N--;

  x.resize(N);
  inf.clear();                    // forget we hit the end of file
  inf.seekg(0, std::ios::beg);    // move to the start of the file
  int i = -1;
  while (inf && ++i < N) {
    if (inf.peek() == '#') inf.ignore(2000,'\n');
    inf >> x[i] >> vi;
    inf.ignore(2000,'\n');
  }
  x.Set();

  Util::getComment(inf,inputf);

  return true;
}
template <typename T, typename M1d>
F1d<T> kramerskronig(const M1d &x, const F1d<T> &f)
{

  F1d<T> g{x.size()};
  #pragma omp parallel for
  for (int i = 0; i < x.size(); i++) {
    double dxinv2 = (i > 0) ? x.d1x(i-1) : 0.0;
    int ip1 = (i < x.size()-1) ? i+1 : i;
    int im1 = (i > 0)    ? i-1 : i;

    //// Somehow below is ~3 times faster...
    //T atSing = (f[ip1]-f[i])*dxInv[i] + (f[i]-f[im1])*dxInv[im1];
    //T sum{0.0};
    //for (int j = 0; j < x.size()-1; j++) {
    //  int jp1 = (j+1 != i) ? j+1 : j;
    //  if (j == i) {
    //    sum += dx*atSing;
    //  } else {
    //    if (jp1 == i) sum += dx*2.0*(f[j]-f[i])/(x[j]-x[i]);
    //    else          sum += dx*((f[j]-f[i])/(x[j]-x[i]) + (f[jp1]-f[i])/(x[jp1]-x[i]));
    //  }
    //}
    //T logVal = (i!=0 && i!=x.size()-1) ? std::log((x.Last()-x[i]) / (x[i]-x.First())) : 0.0;
    //g[i] = (sum/2.0 + logVal)/pi;

    T atSing = 0.5*((f[ip1]-f[i])*x.d1x(i) + (f[i]-f[im1])*dxinv2);
    T sum{0.0};
    for (int j = 0; j < x.size(); j++) {
      if (i!=j) sum += x.dx(j)*(f[j]-f[i])/(x[j]-x[i]);
      else      sum += x.dx(j)*atSing;
    }
    double logVal = (i!=0 && i!=x.size()-1) ? std::log((x.Last()-x[i]) / (x[i]-x.First())) : 0.0;
    g[i] = (sum + f[i]*logVal)/Maths::pi;
  }
  return g;
}
template <typename T>
inline T integralTZ(F1d<T> &f, double dx, int i, int j)
{ // \sum_{k=i}^{j} w_{k}*f_{k}*dx  where  w_{k} (k==i,j)? 1/2 : 1
  T sum{0.0};
  if (j==i)
    return sum;
  sum += 0.5*f[i];
  for (int k = i+1; k < j; k++)
    sum += f[k];
  sum += 0.5*f[j];
  return sum*static_cast<T>(dx);
}
template <typename T>
inline T integralTZ2(F1d<T> &f, double dx, int i, int j)
{ // \sum_{k=i}^{j} w_{k}*f_{k}*dx  where  w_{k} (k==i)? 1/2 : 1
  T sum{0.0};
  //if (j==i)
  //  return sum;
  sum += 0.5*f[i];
  for (int k = i+1; k <= j; k++)
    sum += f[k];
  return sum*static_cast<T>(dx);
}
template <typename T>
template <typename M1d> 
void F1d<T>::cspline(const M1d &X1, const F1d<T> &Y1, const M1d &X2)
{ // Interpolating cubic spline function for irregularly-spaced points.
  // Assumes that X1,X2 entries are monotonically increasing.
  //  Input:
  //      X1      x-coordinates of Y1
  //      Y1      an array of irregular data points (X1.size() entries)
  //      X2      an array of x-coordinates of v(output)
  //  Output:
  //      v       cubic spline sampled at X2 (X2.size() entries)

  assert(X1.First() <= X2.First() && X2.Last() <= X1.Last(), "Can't extrapolate!");
  assert(X1.size() == Y1.size() && X2.size() == N, "Size matters!");

  // compute 1 st derivatives at each point -> YD
  F1d<T> YD{X1.size()};
  getYD_gen(X1, Y1, YD);

  // p1 : left endpoint of interval
  // p2 : resampling position
  // p3 : right endpoint of interval
  // j  : input index for current interval
  double p1{0}, p2{0}, p3{0};
  T A0{0}, A1{0}, A2{0}, A3{0};
  int i, j;
  p3 = X2[0] - 1;                     // force coefficient initialization
  for (i = j = 0; i < X2.size(); i++) {
    // check if in new interval
    p2 = X2[i];
    if (p2 > p3) {
      // find the interval which contains p2
      for (; (j < X1.size() && p2 > X1[j]); j++);
      if (p2 < X1[j]) j--;
      p1 = X1[j];                     // update left endpoint
      p3 = X1[j+1];                   // update right endpoint

      // compute spline coefficients
      double dx = 1.0 / (X1[j+1] - X1[j]);
      T dy = (Y1[j+1] - Y1[j]) * dx;
      A0 = Y1[j];
      A1 = YD[j];
      A2 = dx*   ( 3.0*dy - 2.0*YD[j] - YD[j+1]);
      A3 = dx*dx*(-2.0*dy +     YD[j] + YD[j+1]);
    }

    // use Horner*s rule to calculate cubic polynomial
    double x = p2 - p1;
    v[i] = ((A3*x + A2)*x + A1)*x + A0;
  }
}

template <typename T>
template <typename M1d> 
void F1d<T>::getYD_gen(const M1d &X, const F1d<T> &Y, F1d<T> &YD)
{ // YD <- Computed 1st derivative of data in X,Y (X.size() entries). The not-a-knot boundary condition is used.

  F1d<T> A{X.size()}, B{X.size()}, C{X.size()};

  // init first row data
  double h0 = X[1] - X[0];
  double h1 = X[2] - X[1];
  T r0 = (Y[1] - Y[0]) / h0;
  T r1 = (Y[2] - Y[1]) / h1;
  B[0] = h1*(h0 + h1);
  C[0] = (h0 + h1)*(h0 + h1);
  YD[0] = r0*(3*h0*h1 + 2*h1*h1) + r1*h0*h0;

  // init tridiagonal bands A, B, C, and column Y0
  // YD will later be used to return the derivatives
  int i;
  for (i = 1; i < X.size()-1; i++) {
    h0 = X[i]- X[i-1];
    h1 = X[i+ 1] - X[i];
    r0 = (Y[i]- Y[i-1]) / h0;
    r1 = (Y[i+1]- Y[i]) / h1;
    A[i] = h1;
    B[i] = 2*(h0 + h1);
    C[i] = h0;
    YD[i] = 3*(r0*h1 + r1*h0);
  }

  // last row
  A[i] = (h0 + h1)*(h0 + h1);
  B[i] = h0*(h0 + h1);
  YD[i] = r0*h1*h1 + r1*(3*h0*h1 + 2*h0*h0);

  // solve for the tridiagonal matrix: YD=YD*inv(tridiag matrix)
  tridiag_gen(A,B,C,YD);
}
template <typename T>
void F1d<T>::tridiag_gen(const F1d<T> &A, const F1d<T> &B, const F1d<T> &C, F1d<T> &D)
{ // Gauss Elimination with backsubstitution for general tridiagonal matrix with bands A,B,C and column D.

  int len = A.size();
  F1d<T> F(len);

  // Gauss elimination; forward substitution
  T b = B[0];
  D[0] = D[0] / b;
  int i;
  for (i = 1; i < len; i++) {
    F[i] = C[i-1] / b;
    b = B[i] - A[i]*F[i];
    assert(b != 0.0, "Divide by zero");
    D[i] = (D[i] - D[i-1]*A[i]) / b;
  }
  // backsubstitution
  for(i = len-2; i >= 0; i--) D[i] -= (D[i+1]*F[i+1]);
}

template <typename T>
template <typename Interp> 
inline T F1d<T>::line(const Interp &jd) const
{
  return v[jd.j]+jd.d*(v[jd.j+1]-v[jd.j]);
}


// ============  F2d  ============
// ======  Constructors
template <typename T>
inline F2d<T>::F2d(int Nd_, int N_)
{
  assert0(Nd_ >= 0 && N_ >= 0);
  Nd = Nd_;
  N = N_;

  int offset = sizeof(F1d<T>)*Nd+8;
  mem = operator new(sizeof(T)*Nd*N + offset);
  v = new(mem) F1d<T>[Nd];
  T* begins = reinterpret_cast<T*>(static_cast<char*>(mem)+offset);

  for (int i = 0; i < Nd; i++)
    v[i].init(N,begins+i*N);
}
template <typename T>
inline F2d<T>::F2d(int Nd_, int N_, T vi)
{
  assert0(Nd_ >= 0);
  Nd = Nd_;
  N = N_;

  int offset = sizeof(F1d<T>)*Nd+8;
  mem = operator new(sizeof(T)*Nd*N + offset);
  v = new(mem) F1d<T>[Nd];
  T* begins = reinterpret_cast<T*>(static_cast<char*>(mem)+offset);

  for (int i = 0; i < Nd; i++) {
    v[i].init(N,begins+i*N);
    for (int j = 0; j < N; j++) {
      v[i][j] = vi;
    }
  }
}
template <typename T>
inline F2d<T>::F2d(int Nd_, const F1d<T> &f)
{
  assert0(Nd_ >= 0);
  Nd = Nd_;
  N = f.size();

  int offset = sizeof(F1d<T>)*Nd+8;
  mem = operator new(sizeof(T)*Nd*N + offset);
  v = new(mem) F1d<T>[Nd];
  T* begins = reinterpret_cast<T*>(static_cast<char*>(mem)+offset);

  for (int i = 0; i < Nd; i++) {
    v[i].init(N,begins+i*N);
    for (int j = 0; j < N; j++) {
      v[i][j] = f[j];
    }
  }
}
template <typename T>
inline F2d<T>::F2d(const F2d<T> &F)
{
  Nd = F.Nd;
  N = F.N;

  int offset = sizeof(F1d<T>)*Nd+8;
  mem = operator new(sizeof(T)*Nd*N + offset);
  memcpy(mem, F.mem, sizeof(T)*Nd*N + offset);

  v = new(mem) F1d<T>[Nd];
  T* begins = reinterpret_cast<T*>(static_cast<char*>(mem)+offset);
  for (int i = 0; i < Nd; i++)
    v[i].init(N, begins+i*N);
}

// ====== Operator overloadings {{{
template <typename T>
inline F2d<T> F2d<T>::operator- () const
{
  F2d<T> F{*this};
  for (int i = 0; i < F.Nd; i++)
    for (int j = 0; j < F.N; j++)
      F(i,j) *= static_cast<T>(-1);
  return F;
}
template <typename T>
inline F2d<T>& F2d<T>::operator= (const F2d<T> &F)
{
  if (this == &F)
    return *this;
  if (F.Nd > Nd || F.N > N) {
    Nd = F.Nd;
    N = F.N;

    int offset = sizeof(F1d<T>)*Nd+8;
    operator delete(mem);
    mem = operator new(sizeof(T)*Nd*N + offset);
    memcpy(mem, F.mem, sizeof(T)*Nd*N + offset);

    v = new(mem) F1d<T>[Nd];
    T* begins = reinterpret_cast<T*>(static_cast<char*>(mem)+offset);
    for (int i = 0; i < Nd; i++)
      v[i].init(N, begins+i*N);
  } else {
    Nd = F.Nd;
    N = F.N;
    for (int i = 0; i < Nd; i++)
      memcpy(v[i].begin(), F.v[i].begin(), sizeof(T)*N);
  }

  return *this;
}
template <typename T>
inline F2d<T>& F2d<T>::operator+= (const F2d<T> &F) 
{ 
  assert0((this->Nd == F.Nd) && (this->N == F.N));
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] += F(i,j); 
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator+= (const F1d<T> &f) 
{ 
  assert0(this->N == f.size());
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] += f[j]; 
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator+= (const T &a) 
{ 
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] += a;
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator-= (const F2d<T> &F) 
{ 
  assert0((this->Nd == F.Nd) && (this->N == F.N));
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] -= F(i,j); 
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator-= (const F1d<T> &f) 
{ 
  assert0(this->N == f.size());
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] -= f[j]; 
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator-= (const T &a) 
{ 
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] -= a;
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator*= (const F2d<T> &F) 
{ 
  assert0((this->Nd == F.Nd) && (this->N == F.N));
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] *= F(i,j); 
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator*= (const F1d<T> &f) 
{ 
  assert0(this->N == f.size());
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] *= f[j]; 
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator*= (const T &a) 
{ 
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] *= a;
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator/= (const F2d<T> &F) 
{ 
  assert0((this->Nd == F.Nd) && (this->N == F.N));
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] /= F(i,j); 
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator/= (const F1d<T> &f) 
{ 
  assert0(this->N == f.size());
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] /= f[j]; 
  return *this; 
} 
template <typename T>
inline F2d<T>& F2d<T>::operator/= (const T &a) 
{ 
  assert(a != 0., "Divide by zero(F2d/=a)");
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      (this->v)[i][j] /= a;
  return *this; 
} 

template <typename T>
inline F2d<T> operator+ (const F2d<T> &F1, const F2d<T> &F2)
{
  assert0((F1.Nd == F2.Nd) && (F1.N == F2.N));
  F2d<T> G{F1.Nd, F1.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = F1(i,j) + F2(i,j);
  return G;
} 
template <typename T>
inline F2d<T> operator+ (const F2d<T> &F1, const F1d<T> &f2)
{
  assert0(F1.N == f2.size());
  F2d<T> G{F1.Nd, F1.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = F1(i,j) + f2[j];
  return G;
} 
template <typename T>
inline F2d<T> operator+ (const F1d<T> &f1, const F2d<T> &F2)
{
  return F2 + f1;
} 
template <typename T>
inline F2d<T> operator+ (const F2d<T> &F, const T &a)
{
  F2d<T> G{F.Nd, F.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = F(i,j) + a;
  return G;
} 
template <typename T>
inline F2d<T> operator+ (const T &a, const F2d<T> &F)
{
  return F + a;
} 

template <typename T>
inline F2d<T> operator- (const F2d<T> &F1, const F2d<T> &F2)
{
  assert0((F1.Nd == F2.Nd) && (F1.N == F2.N));
  F2d<T> G{F1.Nd, F1.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = F1(i,j) - F2(i,j);
  return G;
} 
template <typename T>
inline F2d<T> operator- (const F2d<T> &F1, const F1d<T> &f2)
{
  assert0(F1.N == f2.size());
  F2d<T> G{F1.Nd, F1.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = F1(i,j) - f2[j];
  return G;
} 
template <typename T>
inline F2d<T> operator- (const F1d<T> &f1, const F2d<T> &F2)
{
  assert0(f1.size() == F2.N);
  F2d<T> G{F2.Nd, F2.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = f1[j] - F2(i,j);
  return G;
} 
template <typename T>
inline F2d<T> operator- (const F2d<T> &F, const T &a)
{
  return F + (-a);
}
template <typename T>
inline F2d<T> operator- (const T &a, const F2d<T> &F)
{
  F2d<T> G{F.Nd, F.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = a - F(i,j);
  return G;
} 

template <typename T>
inline F2d<T> operator* (const F2d<T> &F1, const F2d<T> &F2)
{
  assert0((F1.Nd == F2.Nd) && (F1.N == F2.N));
  F2d<T> G{F1.Nd, F1.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = F1(i,j) * F2(i,j);
  return G;
} 
template <typename T>
inline F2d<T> operator* (const F2d<T> &F1, const F1d<T> &f2)
{
  assert0(F1.N == f2.size());
  F2d<T> G{F1.Nd, F1.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = F1(i,j) * f2[j];
  return G;
} 
template <typename T>
inline F2d<T> operator* (const F1d<T> &f1, const F2d<T> &F2)
{
  return F2 * f1;
} 
template <typename T>
inline F2d<T> operator* (const F2d<T> &F, const T &a)
{ 
  F2d<T> G{F.Nd, F.N};
  for (int i = 0; i < G.Nd; i++)
    for (int j = 0; j < G.N; j++)
      G(i,j) = F(i,j) * a;
  return G;
} 
template <typename T>
inline F2d<T> operator* (const T &a, const F2d<T> &F)
{
  return F * a;
} 

template <typename T>
inline F2d<T> operator/ (const F2d<T> &F1, const F2d<T> &F2)
{ 
  assert0((F1.Nd == F2.Nd) && (F1.N == F2.N));
  F2d<T> G{F1.Nd, F1.N};
  for (int i = 0; i < G.Nd; i++) {
    for (int j = 0; j < G.N; j++) {
      assert(F2(i,j) != 0., "Divide by zero");
      G(i,j) = F1(i,j) / F2(i,j);
    }
  }
  return G;
} 
template <typename T>
inline F2d<T> operator/ (const F2d<T> &F1, const F1d<T> &f2)
{ 
  assert0(F1.N == f2.size());
  F2d<T> G{F1.Nd, F1.N};
  for (int i = 0; i < G.Nd; i++) {
    for (int j = 0; j < G.N; j++) {
      assert(f2[i] != 0., "Divide by zero");
      G(i,j) = F1(i,j) / f2[j];
    }
  }
  return G;
} 
template <typename T>
inline F2d<T> operator/ (const F1d<T> &f1, const F2d<T> &F2)
{ 
  assert0(f1.size() == F2.N);
  F2d<T> G{F2.Nd, F2.N};
  for (int i = 0; i < G.Nd; i++) {
    for (int j = 0; j < G.N; j++) {
      assert(F2(i,j) != 0., "Divide by zero");
      G(i,j) = f1[j] / F2(i,j);
    }
  }
  return G;
} 
template <typename T>
inline F2d<T> operator/ (const F2d<T> &F, const T &a)
{
  assert(a != 0., "Divide by zero");
  return F * (1./a);
} 
template <typename T>
inline F2d<T> operator/ (const T &a, const F2d<T> &F)
{
  F2d<T> G{F.Nd, F.N};
  for (int i = 0; i < G.Nd; i++) {
    for (int j = 0; j < G.N; j++) {
      assert(F(i,j) != 0., "Divide by zero");
      G(i,j) = a / F(i,j);
    }
  }
  return G;
}
// }}}

template <typename T>
inline void F2d<T>::resize(int Nd_, int N_)
{
  assert0(Nd_ >= 0 && N_ >= 0);
  if (Nd_ > Nd || N_>N) {
    operator delete(mem);

    Nd = Nd_;
    N = N_;

    int offset = sizeof(F1d<T>)*Nd+8;
    mem = operator new(sizeof(T)*Nd*N + offset);
    v = new(mem) F1d<T>[Nd];

    T* begins = reinterpret_cast<T*>(static_cast<char*>(mem)+offset);
    for (int i = 0; i < Nd; i++)
      v[i].init(N, begins+i*N);
  } else {
    Nd = Nd_;
    N = N_;
  }
}
template <typename T>
inline void F2d<T>::Set(int Nd_, int N_, T vi)
{
  resize(Nd_, N_);
  for (int i = 0; i < Nd; i++)
    for (int j = 0; j < N; j++)
      v[i][j] = vi;
}
template <typename T>
inline double F2d<T>::Rms()
{
  double sum{0.0};
  for (int i = 0; i < Nd; i++)
    sum += v[i].Rms();
  return sum/Nd;
}
template <typename T>
double normalizedDiff(const F2d<T> &F1, const F2d<T> &F2)
{
  assert0(F1.Nd == F2.Nd);
  double ndif{0.};
  for (int i = 0; i < F1.Nd; i++)
    ndif += normalizedDiff(F1[i],F2[i]);
  return ndif/F1.Nd;
}
template <typename T>
template <typename M1d> 
void F2d<T>::Print(const M1d &x, const std::string &outputf, const std::string &comment)
{
  assert0(x.size() == N);
  std::clog << "  => Function printed: " << outputf << "\n";

  std::ofstream outf{outputf};
  if (comment != "")
    outf << "# " << comment << '\n';
  outf << std::fixed << std::setprecision(25);
  for (int i = 0; i < N; i++) {
    outf << x[i] << '\t';
    for (int j = 0; j < Nd; j++) {
      outf << v[j][i] << '\t';
    }
    outf << '\n';
  }
}
template <typename T>
bool F2d<T>::Read(const std::string &inputf, const std::string &descrip)
{
  std::clog << "Reading " << descrip << " from " << inputf << ", with " << Nd <<"x" << N << " entries.\n";

  std::ifstream inf{inputf};
  if (!inf.good()) {
    std::clog << "Missing input file! " << "No " << inputf << "\n";
    return false;
  }

  double x;
  int i = -1;
  while (inf && ++i < N) {
    if (inf.peek() == '#') inf.ignore(2000,'\n');
    inf >> x;
    for (int j = 0; j < Nd; j++)
      inf >> v[j][i];
    inf.ignore(2000,'\n');
  }
  if (i!=N) {
    std::clog << "Reading failed: " << inputf << " i(=" << i << ") != N (" << N << ")\n";
    return false;
  }
  Util::getComment(inf,inputf);

  return true;
}
template <typename T>
bool F2d<T>::Readr(const std::string &inputf)
{
  std::clog << "Reading " << inputf << ", with " << Nd <<"x" << N << " entries in the row major fashion.\n";

  std::ifstream inf{inputf};
  if (!inf.good()) {
    std::clog << "Missing input file! " << "No " << inputf << "\n";
    return false;
  }

  int i = -1;
  while (inf && ++i < N) {
    if (inf.peek() == '#') inf.ignore(2000,'\n');
    for (int j = 0; j < Nd; j++)
      inf >> v[i][j];
    inf.ignore(2000,'\n');
  }
  if (i!=N) {
    std::clog << "Reading failed: " << inputf << " i(=" << i << ") != N (" << N << ")\n";
    return false;
  }
  Util::getComment(inf,inputf);

  return true;
}
template <typename T>
template <typename M1d> 
bool F2d<T>::Readx(M1d &x, const std::string &inputf, const std::string &descrip)
{
  std::clog << "Reading " << descrip << " with mesh from " << inputf << "\n";
  if (Nd <= 0) {
    std::clog << "Readx: Nd must be set before.\n";
    return false;
  }

  std::ifstream inf{inputf};
  if (!inf.good()) {
    std::clog << "Missing input file! " << "No " << inputf << "\n";
    return false;
  }

  this->resize(Nd,0);
  double xi;
  T vij;
  while (inf && ++N) {
    if (inf.peek() == '#') inf.ignore(2000,'\n');
    inf >> xi;
    for (int i = 0; i < Nd; i++) {
      inf >> vij;
      v[i].push_back(vij);
    }
    inf.ignore(2000,'\n');
  }
  N--;

  x.resize(N);
  inf.clear();               // forget we hit the end of file
  inf.seekg(0, std::ios::beg);    // move to the start of the file
  int i = -1;
  while (inf && ++i < N) {
    if (inf.peek() == '#') inf.ignore(2000,'\n');
    inf >> x[i];
    for (int j = 0; j < Nd; j++) {
      inf >> vij;
    }
    inf.ignore(2000,'\n');
  }
  x.Set();

  Util::getComment(inf,inputf);

  return true;
}
template <typename T, typename M1d>
F2d<T> kramerskronig(const M1d &x, const F2d<T> &F)
{
  F2d<T> G{F.Nd, F.N};
  for (int i = 0; i < G.Nd; i++)
    G[i] = kramerskronig(x, F[i]);
  return G;
}


#endif /* end of include guard: FUNCTION_H */
