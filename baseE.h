// 17 Oct 2018 16:07:38 -0400
#ifndef BASEE_H
#define BASEE_H

#include <cmath>
#include <iostream>

class baseE
{ // Real number represented as r e^t
  private:
    double r, t;
  public:
    baseE(double r_ = 0, double t_ = 0) : r(r_), t(t_) { }
    baseE(const baseE& z) : r(z.r), t(z.t) { }
    virtual ~baseE() { }

    //friend baseE operator- (const baseE& z) { return baseE(-z.r,z.t);};
    baseE operator- () const { return baseE(-r, t); }

    baseE& operator=  (const baseE& z);
    baseE& operator+= (const baseE& z);
    baseE& operator-= (const baseE& z);
    baseE& operator*= (const baseE& z);
    baseE& operator/= (const baseE& z);

    friend bool operator== (const baseE& x, const baseE& y) { return x.r == y.r && x.t == y.t; }
    friend bool operator!= (const baseE& x, const baseE& y) { return !(x == y); }
    friend bool operator>  (const baseE& x, const baseE& y);
    friend bool operator<  (const baseE& x, const baseE& y);
    friend bool operator>= (const baseE& x, const baseE& y) { return !operator<(x, y); }
    friend bool operator<= (const baseE& x, const baseE& y) { return !operator>(x, y); }

    friend bool operator== (const baseE& x, double y) { return x.r*exp(x.t) == y; }
    friend bool operator== (double x, const baseE& y) { return y == x; }
    friend bool operator!= (const baseE& x, double y) { return !(x == y); }
    friend bool operator!= (double x, const baseE& y) { return !(x == y); }

    friend baseE operator+ (const baseE& x, const baseE& y);
    friend baseE operator- (const baseE& x, const baseE& y);
    friend baseE operator* (const baseE& x, const baseE& y) { return baseE(x.r*y.r, x.t+y.t); }
    friend baseE operator/ (const baseE& x, const baseE& y) { return baseE(x.r/y.r, x.t-y.t); }

    friend baseE operator+ (const baseE& x, double y) { return x + baseE(y); }
    friend baseE operator+ (double x, const baseE& y) { return baseE(x) + y; }
    friend baseE operator- (const baseE& x, double y) { return x - baseE(y); }
    friend baseE operator- (double x, const baseE& y) { return baseE(x) - y; }
    friend baseE operator* (const baseE& x, double y) { return baseE(x.r*y, x.t); }
    friend baseE operator* (double x, const baseE& y) { return baseE(x*y.r, y.t); }
    friend baseE operator/ (const baseE& x, double y) { return baseE(x.r/y, x.t); }
    friend baseE operator/ (double x, const baseE& y) { return baseE(x/y.r,-y.t); }
    friend double div(const baseE& x, const baseE& y) { return x.r/y.r*exp(x.t-y.t); }
  
    friend std::ostream& operator<< (std::ostream& oS, const baseE& z) { oS << z.r << " " << z.t << " "; return oS; }

    baseE& balance();
    double& mantissa() { return r; }
    double& exponent() { return t; }
    double mantissa() const { return r; }
    double exponent() const { return t; }
    double decimal() const { return r*exp(t); }
    //double exp_dbl() const { return log(abs(r))+t; }

    friend baseE sqrt(const baseE& z) { return baseE(sqrt(z.r), z.t/2.0); }
    friend baseE abs(const baseE& z)  { return baseE(fabs(z.r), z.t); }
    friend bool isnan(const baseE& z) { return std::isnan(z.r) || std::isnan(z.t); }
    friend double log(const baseE& z) { return log(z.r) + z.t; }
};

inline baseE& baseE::operator=  (const baseE& z)
{
  r = z.r;
  t = z.t;
  return *this;
}

inline baseE& baseE::operator+= (const baseE& z)
{
  if (z.t<t)
    r += z.r*exp(z.t-t);
  else {
    r = z.r+r*exp(t-z.t);
    t=z.t;
  }
  return *this;
}

inline baseE& baseE::operator-= (const baseE& z)
{
  if (z.t<t)
    r -= z.r*exp(z.t-t);
  else{
    r = -z.r+r*exp(t-z.t);
    t=z.t;
  }
  return *this;
}

inline baseE& baseE::operator*= (const baseE& z)
{
  r *= z.r;
  t += z.t;
  return *this;
}

inline baseE& baseE::operator/= (const baseE& z)
{
  r /= z.r;
  t -= z.t;
  return *this;
}

inline bool operator> (const baseE& x, const baseE& y)
{
  return (x.t > y.t) ? (x.r > y.r * exp(y.t-x.t)) : (x.r * exp(x.t-y.t) > y.r);
}
inline bool operator< (const baseE& x, const baseE& y)
{
  return (x.t > y.t) ? (x.r < y.r * exp(y.t-x.t)) : (x.r * exp(x.t-y.t) < y.r);
}

inline baseE operator+ (const baseE& x, const baseE& y)
{
  return (x.t>y.t) ? baseE(x.r+y.r*exp(y.t-x.t), x.t) : baseE(y.r+x.r*exp(x.t-y.t), y.t);
}

inline baseE operator- (const baseE& x, const baseE& y)
{
  return (x.t>y.t) ? baseE(x.r-y.r*exp(y.t-x.t), x.t) : baseE(-y.r+x.r*exp(x.t-y.t), y.t);
}

inline baseE& baseE::balance()
{
  if (r==0) {
    t=0;
    return *this;
  }
  t += log(fabs(r));
  r = (r>0) ? 1 : -1;
  return *this;
}

//baseE balance(const baseE& x)
//{
//  if (x.r == 0) return baseE(0.0, 0.0);
//  return baseE( (x.r > 0) ? 1. : -1., x.t + log(fabs(x.r)) );
//}
#endif /* end of include guard: BASEE_H */
