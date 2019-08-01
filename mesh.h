#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "maths.h"
#include "assert.h"
#include "zeroin.h"
#include "util.h"

class M1d;
class AddLogTan;
class Interp;

class M1d
{
  protected:
    int N;
    int dN;
    double *x;    // x = x_i, i \in [0,N-1]
    double *Dx;   // \Delta x = \frac{x_{i+1} - x_{i-1}}{2}, i \in [1,N-2]
                  //    (Dx[0] = \frac{x_1 - x_0}{2}, Dx[N-1] = \frac{x_{N-1} - x_{N-2}}{2}) for trapaziod )
    double *D1x;  // \frac{1}{\Delta x} = \frac{1}{x_{i+1} - x_{i}}, x \in [0,N-2]
                  //    (D1x[N-1] = 0 not to extrapolate)
  
  public:
    M1d() : N{0}, dN{0}, x{nullptr}, Dx{nullptr}, D1x{nullptr} { }
    M1d(int N_);
    M1d(const M1d &m);
    virtual ~M1d();

    M1d& operator= (const M1d &m);

    double& operator[] (const int i) { assert(0<=i && i<N, "i="<<i<<" N="<<N); return x[i]; }
    const double& operator[] (const int i) const { assert(0<=i && i<N, "i="<<i<<" N="<<N); return x[i]; }
   
    double* begin() const { return x; }
    double* end()   const { return x+N; }

    double& First() { return x[0]; }
    double& Last()  { return x[N-1]; }
    const double& First() const { return x[0]; }
    const double& Last()  const { return x[N-1]; }

    double dx(const int i) const { assert(0<=i && i<N, "i="<<i<<" N="<<N); return Dx[i]; }
    double d1x(const int i) const { assert(0<=i && i<N, "i="<<i<<" N="<<N); return D1x[i]; }

    int size() const { return N; }
    void resize(int N_ = 0);
    void Set();

    void makeEqDist(int N_, double a, double b);
    void makeTan   (int N_, double dx, double x1); // [-x1, x1]
    void makeLogTan(int N_, double x0, double x1, double x2, double alpha); // [x0, x2)

    void makeEqLogTan    (int N_, double x0, double x1, double x2, double alpha); // [-x2 ~ x2]
    void makePosEqLogTan (int N_, double x0, double x1, double x2, double alpha); // (0.0 ~ x2)

    void makeEqLogTan0   (int N_, double x0, double x1, double x2, double alpha); // (-x2 ~ 0 ~ x2) 
    void makePosEqLogTan0(int N_, double x0, double x1, double x2, double alpha); // [0.0 ~ x2)

    int fInd(double a, int& i) const;
    int finD(double a, int& i) const;
    int fInd(double a) const;
    int finD(double a) const;

    int whereIs(double a, int &i) const;
    int whereIs(double a) const;
    int whereIsD(double a, int &i) const;
    int whereIsD(double a) const;

    using locateFcn = int (M1d::*)(double,int&) const;
    Interp getInterp(const double a, int &i, locateFcn locate) const;

    bool Read(const std::string &inputf);
    void Print(const std::string &outputf, const std::string &comment = "");

  private:
    int bisection(double a, int &jl, int &ju) const;
};

class AddLogTan
{ // functor class for adding log mesh and tan mesh using zeroin.
  private:
    double dwt, xt;
  public:
    static const double solMin, solMax, solTol;
    
    AddLogTan(double dwt_, double xt_) : dwt{dwt_}, xt{xt_} { }
    virtual ~AddLogTan() { }
    double operator() (double u);
};

class Interp {
  public:
    int j;
    double d;
    Interp() : j{0}, d{0.0} { }
    Interp(int j_, double d_) : j{j_}, d{d_} { }
    virtual ~Interp() { }
};


// ============  M1d  ============
inline M1d::M1d(int N_)
{
  assert0(N_ >= 0);
  N = N_;
  dN = Maths::max(10,static_cast<int>(std::pow(N,0.25)));
  x = new double[N];
  Dx = new double[N];
  D1x = new double[N];
}
inline M1d::M1d(const M1d &m)
{
  N = m.N;
  dN = Maths::max(10,static_cast<int>(std::pow(N,0.25)));
  x = new double[N];
  Dx = new double[N];
  D1x = new double[N];
  std::copy(m.x, m.x+m.N, x);
  std::copy(m.Dx, m.Dx+m.N, Dx);
  std::copy(m.D1x, m.D1x+m.N, D1x);
}
inline M1d::~M1d()
{
  delete[] x;
  delete[] Dx;
  delete[] D1x;
  x = nullptr;
  Dx = nullptr;
  D1x = nullptr;
}
inline M1d& M1d::operator= (const M1d &m)
{
  if (this == &m)
    return *this;
  this->N = m.N;
  this->dN = m.dN;
  if (this->x) delete[] this->x;
  if (this->Dx) delete[] this->Dx;
  if (this->D1x) delete[] this->D1x;
  this->x = new double[this->N];
  this->Dx = new double[this->N];
  this->D1x = new double[this->N];
  std::copy(m.x, m.x+m.N, this->x);
  std::copy(m.Dx, m.Dx+m.N, this->Dx);
  std::copy(m.D1x, m.D1x+m.N, this->D1x);
  return *this;
}
inline void M1d::resize(int N_)
{
  assert0(N_ >= 0);
  if (N_ > N) {
    if (x) delete[] x;
    if (Dx) delete[] Dx;
    if (D1x) delete[] D1x;
    x = new double[N_];
    Dx = new double[N_];
    D1x = new double[N_];
  }
  N = N_;
  dN = Maths::max(10,static_cast<int>(std::pow(N,0.25)));
}
inline void M1d::Set()
{
  Dx[0] = 0.5*(x[1]-x[0]);
  D1x[0] = 1./(x[1]-x[0]);
  for (int i = 1; i < N-1; i++) {
    Dx[i] = 0.5*(x[i+1]-x[i-1]);
    D1x[i] = 1./(x[i+1]-x[i]);
  }
  Dx[N-1] = 0.5*(x[N-1]-x[N-2]);
  D1x[N-1] = 0.0;
}
inline void M1d::makeEqDist(int N_, double a, double b)
{
  resize(N_);
  for (int i = 0; i < N; i++) {
    x[i] = a + (b-a)/(N-1)*i;
  }
  Set();
}
inline void M1d::makeTan(int N_, double dx, double x1)
{ // x \in [-x1, x1]
  // Note: 0.0 is (not) included if N is (even) odd
  resize(N_);
  double piN = Maths::pi/N;
  double dxx1 = dx/x1;
  auto f = [&](double d) { return std::tan(d)*std::tan(piN-2.0*d/N) - dxx1; };
  double d = zeroin(1e-7, Maths::pi/2.0 - 1e-7, f, 1e-15);
  double w = x1*std::tan(d);
  double a = Maths::pi/2.0-d;

  x[0] = -x1;
  for (int i = 1; i < N-1; i++)
    x[i] = w*std::tan(a*(2.0*i/(N-1) - 1.0));
  x[N-1] = x1;
  Set();
}

inline void M1d::makeLogTan(int N_, double x0, double x1, double x2, double alpha)
{ // x \in [x0, x2), x0 > 0
  //  1) x0 ~ x1 : logarithm mesh,
  //  2) x1 ~ x2 : tangent mesh
  // Note: x1/x0 ~ 1e5,  x2/x1 ~ 10, alpha ~ 0.2 seems enough
  assert(x0 < x1 && x1 < x2, "x0,x1,x2="<<x0<<","<<x1<<","<<x2);
  assert(alpha >= 0, "alpha="<<alpha);
  resize(N_);
  double Dlogx = std::log(x1/x0);
  double eta = Dlogx/(x2/x1-1);
  double N1_min = (1+eta*N)/(1+eta)+0.5;
  int N1 = static_cast<int>((1+alpha)*N1_min);
  if (N1 > N-2)
    N1 = N-2;
  
  int N2 = N-N1;
  double xt = x2/x1;
  double dwt = N2*Dlogx/(N1-1);

  AddLogTan f(dwt,xt);
  double ut = zeroin(AddLogTan::solMin,AddLogTan::solMax,f,AddLogTan::solTol);
  
  double a = std::atan(std::tan(ut)/xt);
  double b = dwt*std::sin(a)*std::cos(a);
  double w = x1/std::tan(a);

  for (int i = 0; i < N1; i++)
    x[i] = x0*std::exp(i*Dlogx/(N1-1));
  for (int i = 0; i < N2; i++)
    x[N1 + i] = w*std::tan(a+(i+1)*b/N2);
  Set();
}
inline void M1d::makePosEqLogTan0(int N_, double x0, double x1, double x2, double alpha)
{ // x \in [0, x2)
  //  1) 0. ~ x0 : equi-distance mesh,
  //  2) x0 ~ x1 : logarithm mesh,
  //  3) x1 ~ x2 : tangent mesh
  assert(x0 < x1 && x1 < x2, "x0,x1,x2="<<x0<<","<<x1<<","<<x2);
  assert(alpha >= 0, "alpha="<<alpha);
  resize(N_);
  double Dlogx = std::log(x1/x0);
  double eta = Dlogx/(x2/x1-1);
  double N1_min = (1+eta*(N+1/Dlogx))/(1+eta*(1+1/Dlogx))+0.5;
  int N1 = static_cast<int>((1+alpha)*N1_min);
  int N0 = static_cast<int>((N1-1)/Dlogx);
  if (N1 > N-4) {
    N1 = N-4;
    N0 = 2;
  } else if (N0 > N-N1-2) {
    N0 = N-N1-2;
  }
  
  int N2 = N-N1-N0;
  double xt = x2/x1;
  double dwt = N2*Dlogx/(N1-1);

  AddLogTan f(dwt,xt);
  double ut = zeroin(AddLogTan::solMin,AddLogTan::solMax,f,AddLogTan::solTol);
  
  double a = std::atan(std::tan(ut)/xt);
  double b = dwt*std::sin(a)*std::cos(a);
  double w = x1/std::tan(a);

  for (int i = 0; i < N0; i++)
    x[i] = i*x0/N0;
  for (int i = 0; i < N1; i++)
    x[N0+i] = x0*std::exp(i*Dlogx/(N1-1));
  for (int i = 0; i < N2; i++)
    x[N0+N1+i] = w*std::tan(a+(i+1)*b/N2);
  Set();
}
inline void M1d::makePosEqLogTan(int N_, double x0, double x1, double x2, double alpha)
{ // x \in (0, x2)
  //  1) 0. ~ x0 : equi-distance mesh,
  //  2) x0 ~ x1 : logarithm mesh,
  //  3) x1 ~ x2 : tangent mesh
  assert(x0 < x1 && x1 < x2, "x0,x1,x2="<<x0<<","<<x1<<","<<x2);
  assert(alpha >= 0, "alpha="<<alpha);
  resize(N_);
  double Dlogx = std::log(x1/x0);
  double eta = Dlogx/(x2/x1-1);
  double N1_min = (1+eta*(N+1/Dlogx))/(1+eta*(1+1/Dlogx))+0.5;
  int N1 = static_cast<int>((1+alpha)*N1_min);
  int N0 = static_cast<int>((N1-1)/Dlogx);
  if (N1 > N-4) {
    N1 = N-4;
    N0 = 2;
  } else if (N0 > N-N1-2) {
    N0 = N-N1-2;
  }
  
  int N2 = N-N1-N0;
  double xt = x2/x1;
  double dwt = N2*Dlogx/(N1-1);

  AddLogTan f(dwt,xt);
  double ut = zeroin(AddLogTan::solMin,AddLogTan::solMax,f,AddLogTan::solTol);
  
  double a = std::atan(std::tan(ut)/xt);
  double b = dwt*std::sin(a)*std::cos(a);
  double w = x1/std::tan(a);

  for (int i = 0; i < N0; i++)
  k x[i] = (2*i+1)*x0/(2*N0+1);
  for (int i = 0; i < N1; i++)
    x[N0+i] = x0*std::exp(i*Dlogx/(N1-1));
  for (int i = 0; i < N2; i++)
    x[N0+N1+i] = w*std::tan(a+(i+1)*b/N2);
  Set();
}
inline void M1d::makeEqLogTan0(int N_, double x0, double x1, double x2, double alpha)
{ // x \in (-x2, x2)
  //  1) -x0 ~  x0 : equi-distance mesh,
  //  2) |x0 ~ x1| : logarithm mesh,
  //  3) |x1 ~ x2| : tangent mesh
  // Note: 0.0 is included.
  int n = N_/2+1;
  makePosEqLogTan0(n,x0,x1,x2,alpha);

  double *tmp = new double[n];
  std::copy(x,x+n,tmp);

  resize(N_);
  for (int i = 0; i < n-1; i++)
    x[i] = -tmp[n-1-i];
  for (int i = 0; i < n; i++)
    x[n-1+i] = tmp[i];
  delete[] tmp;
  Set();
}
inline void M1d::makeEqLogTan(int N_, double x0, double x1, double x2, double alpha)
{ // x \in [-x2, x2]
  //  1) -x0 ~  x0 : equi-distance mesh,
  //  2) |x0 ~ x1| : logarithm mesh,
  //  3) |x1 ~ x2| : tangent mesh
  // Note: 0.0 is not included.
  int n = N_>>1;
  makePosEqLogTan(n,x0,x1,x2,alpha);

  double *tmp = new double[n];
  std::copy(x,x+n,tmp);
  tmp[n-1] = x2;

  resize(n<<1);
  for (int i = 0; i < n; i++) {
    x[i] = -tmp[n-1-i];
    x[n+i] = tmp[i];
  }
  delete[] tmp;
  Set();
}
inline int M1d::whereIs(double a, int &i) const
{ // returns i where, x[i] <= a < x[i+1], i = [0,N-1). i is found by increasing it
  if (a < x[i+1]) return i;
  while (i+1 < N-1 && x[i+1] <= a)
    i++;
  return i;
}
inline int M1d::whereIs(double a) const
{
  int i = 0;
  return whereIs(a, i);
}
inline int M1d::whereIsD(double a, int &i) const
{ // returns i where, x[i] <= a < x[i+1], i = [0,N-1). i is found by decreasing it
  if (a > x[i]) return i;
  while (i > 0 && x[i] > a)
    i--;
  return i;
}
inline int M1d::whereIsD(double a) const
{
  int i = N-1;
  return whereIsD(a, i);
}
inline int M1d::bisection(double a, int &jl, int &ju) const
{
  int j = 0;
  while (ju-jl > 1) {
    j = (ju+jl)>>1;
    if (x[j] <= a) 
      jl=j;
    else
      ju=j;
  }
  return jl;
}
inline int M1d::fInd(double a, int& i) const
{
  if (a < x[i+1])
    return i;
  if (i >= N-2)
    return i;

  i++;
  if (a < x[i+1])
    return i;
  if (i >= N-2) 
    return i;

  i++;
  int iu = N-1;
  if (i+dN > N){
    return bisection(a, i, iu);
  }
  if (a < x[i+dN-1]) {
    iu = i+dN-1;
    return bisection(a, i, iu);
  }

  i += dN-2;
  if (i >= N-2) 
    return bisection(a, i, iu);

  i++;
  return bisection(a, i, iu);
}
inline int M1d::fInd(double a) const
{
  int i = 0;
  return fInd(a,i);
}
inline int M1d::finD(double a, int& i) const
{
  if (a >= x[i])
    return i;
  if (i<=0)
    return i;

  i--;
  if (a >= x[i])
    return i;
  if (i-1 <= 0)
    return (i = 0);
  if (a >= x[i-1])
    return --i;

  int il = 0;
  if (i-dN < 0)
    return i = bisection(a, il, i); 
  if (a >= x[i-dN+1]) {
    il = i-dN+1;
    return i = bisection(a, il, --i); 
  }

  i -= dN-2;
  return i = bisection(a, il, i);
}
inline int M1d::finD(double a) const
{
  int i = N-2;
  return finD(a,i);
}
inline Interp M1d::getInterp(const double a, int &i, locateFcn locate) const
{
  int j = (this->*locate)(a, i);
  double d = (a-x[j])*D1x[j];
  return Interp{j,d};
}
bool M1d::Read(const std::string &inputf)
{
  std::ifstream inf{inputf};
  if (!inf.good()) {
    std::clog << "Missing input file! " << "No " << inputf << "\n";
    return false;
  }

  int n = std::count(std::istreambuf_iterator<char>(inf), std::istreambuf_iterator<char>(), '\n');

  inf.clear();                    // forget we hit the end of file
  inf.seekg(0, std::ios::beg);    // move to the start of the file
  if (inf.peek() == '#')
    n--;

  resize(n);

  int i = -1;
  while (inf && ++i < N) {
    if (inf.peek() == '#') inf.ignore(2000,'\n');
    inf >> x[i];
    inf.ignore(2000,'\n');
  }
  if (i!=N) {
    std::clog << "Reading failed: " << inputf << " i(=" << i << ") != N (" << N << ")\n";
    return false;
  }
  Util::getComment(inf,inputf);

  return true;
}
void M1d::Print(const std::string &outputf, const std::string &comment)
{
  std::ofstream outf{outputf};
  if (comment != "")
    outf << "# " << comment << '\n';
  
  outf << "# i\tx[i]\t\t\t Dx[i]\t\t\tD1x[i]\n";

  outf << std::fixed << std::setprecision(15);
  for (int i = 0; i < N; i++)
    outf << i << "\t" << x[i] << std::setw(24) << Dx[i] << std::setw(24) << D1x[i] << '\n';
}

// ============  AddLogTan  ============
const double AddLogTan::solMin = 1e-7;
const double AddLogTan::solMax = Maths::pi/2 - 1e-7;
const double AddLogTan::solTol = 1e-10;
inline double AddLogTan::operator() (double u)
{
  double tg = std::tan(u);
  return u-std::atan(tg/xt)-dwt*xt*tg/(Maths::sqr(xt)+Maths::sqr(tg));
}


#endif /* end of include guard: MESH_H */
