#ifndef MATHS_H
#define MATHS_H

#include <cmath>
#include <limits>

namespace Maths
{
  // constants
  extern const double pi         = 3.14159265358979323846264338327950288419716939937511;
//extern const double e          = 2.71828182845904523536028747135266249775724709369995;
  extern const double Inf        = std::numeric_limits<double>::infinity();
  extern const double NaN        = std::numeric_limits<double>::quiet_NaN();
  extern const double absEpsilon = std::numeric_limits<double>::epsilon();
  extern const double doubleMin  = std::numeric_limits<double>::min();
  extern const double doubleMax  = std::numeric_limits<double>::max();
  extern const double relEpsilon = 1.0e-11;

  bool isEqual(double a, double b, double rE = relEpsilon, double aE = absEpsilon);

  template <typename T> int sgn(T a);
  template <typename T> T sqr(T a);
  template <typename T> T sign(const T &a, const T &b);
  template <typename T> void swap(T &a, T &b);
  template <typename T> inline const T& max(const T &a, const T &b);
  template <typename T> inline const T& min(const T &a, const T &b);
  template <typename T> inline T abs(T a);
  template <typename T> inline T delta(const int i, const int j);

  int binomial(int n, int m);
  class QuadEq;
}

inline bool Maths::isEqual(double a, double b, double rE, double aE)
{
  double diff = abs(a - b);
  return diff <= aE ? 1 : diff <= rE*Maths::max(abs(a),abs(b));
}

template <typename T> inline int Maths::sgn(T a)
{
  return (T{0} < a) - (a < T{0});
}

template <typename T> inline T Maths::sqr(T a)
{
  return a * a;
}

template <typename T> inline T Maths::sign(const T &a, const T &b)
{
  return (b >= 0) ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

template <typename T> inline void Maths::swap(T &a, T &b)
{
  T c = a; a = b; b = c;
}

template <typename T> inline const T& Maths::max(const T &a, const T &b)
{
  return b > a ? b : a;
}

template <typename T> inline const T& Maths::min(const T &a, const T &b)
{
  return b < a ? b : a;
}

template <typename T> inline T Maths::abs(T a)
{
  return a < 0 ? -a : a;
} 

template <typename T> inline T Maths::delta(const int i, const int j)
{
  return (i==j) ? T(1) : T(0);
}

int Maths::binomial(int n, int m)
{
  double A = 1.0;
  for (int i = n; i >= n-m+1; i--) A *= i;
  double B = 1.0;
  for (int i = 2; i <= m; i++)     B *= i;

  return static_cast<int>(A/B);
}


class Maths::QuadEq
{ // Quadratic equation : Q(x) = C_0 (x-x0)^2 + 2 C_1 (x-x0) + C_2
  private:
    const double C[3];
    const double x0;

  public:
    QuadEq(double a, double b, double c, double x0_) : C{a,b,c}, x0{x0_} { }
    virtual ~QuadEq() { }
    double operator()(const double x);
    int solve(double &x1, double &x2);
};

int Maths::QuadEq::solve(double &x1, double &x2)
{ // returns number of solutions of Q(x) = 0, and solutions x1, x2
  double Discriminant = sqr(C[1]) - C[0]*C[2];
  if (Discriminant < 0) {
    x1 = x2 = NaN;
    return 0;
  }
  double q = -(C[1]+sgn(C[1])*Discriminant);
  x1 = q/C[0] + x0;
  x2 = C[2]/q + x0;
  if (Discriminant == 0)
    return 1;
  return 2;
}

double Maths::QuadEq::operator()(const double x)
{ // returns Q(x)
  double xx = x - x0;
  return xx*(xx*C[0] + 2.0*C[1]) + C[2];
}

#endif /* end of include guard: MATHS_H */
