#ifndef COMPLEX_H
#define COMPLEX_H
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

class compdb
{ // Custom complex number.
  private:
    double re, im;
  public:
    compdb(double r = 0, double i = 0) : re(r), im(i) { }
    compdb(const compdb &z) : re(z.re), im(z.im) { }
    ~compdb() { } // No virtual here: it adds 8 bytes so that sizeof(compdb) = 24.

    compdb operator- () const;

    compdb& operator=  (const compdb &z);
    compdb& operator=  (const double r);
    compdb& operator+= (const compdb &z);
    compdb& operator+= (const double r);
    compdb& operator-= (const compdb &z);
    compdb& operator-= (const double r);
    compdb& operator*= (const compdb &z);
    compdb& operator*= (const double r);
    compdb& operator/= (const compdb &z);
    compdb& operator/= (const double r);

    friend bool operator== (const compdb &x, const compdb &y);
    friend bool operator!= (const compdb &x, const compdb &y);

    friend bool operator== (const compdb &x, double y);
    friend bool operator== (double x, const compdb &y);
    friend bool operator!= (const compdb &x, double y);
    friend bool operator!= (double x, const compdb &y);
    friend compdb operator+ (const compdb &x, const compdb &y);
    friend compdb operator- (const compdb &x, const compdb &y);
    friend compdb operator* (const compdb &x, const compdb &y);
    friend compdb operator/ (const compdb &x, const compdb &y);
    friend compdb operator+ (const compdb &x, double y);
    friend compdb operator+ (double x, const compdb &y);
    friend compdb operator- (const compdb &x, double y);
    friend compdb operator- (double x, const compdb &y);
    friend compdb operator* (const compdb &x, double y);
    friend compdb operator* (double x, const compdb &y);
    friend compdb operator/ (const compdb &x, double y);
    friend compdb operator/ (double x, const compdb &y);

    friend std::ostream& operator<<(std::ostream &oS, const compdb &z);
    friend std::istream& operator>>(std::istream &iS, compdb &z);

    double& real() { return re; }
    double& imag() { return im; }
    double real() const { return re; }
    double imag() const { return im; }
    compdb conj() const { return compdb(re,-im); }


    friend double real(const compdb &z) { return z.re; }
    friend double imag(const compdb &z) { return z.im; }
    friend double& real(compdb &z) { return z.re; }
    friend double& imag(compdb &z) { return z.im; }
    friend compdb conj(const compdb &z) { return compdb(z.re, -z.im); }
    //friend compdb conj(const double &r) { return compdb{r,0}; }
    friend compdb polar(double r, double phi);

    friend compdb sqrt(const compdb &z);
    friend compdb sqrt_(const compdb &z);
    friend compdb sqrt(const compdb &z1, const compdb &z2);
    friend double abs(const compdb &z);
    friend double norm(const compdb &z) { return z.re*z.re + z.im*z.im; }
    friend bool isnan(const compdb &z) { return std::isnan(z.re) || std::isnan(z.im); }
    friend double arg(const compdb &z) { return atan2(z.im, z.re); }
    friend compdb exp(const compdb &z);
    friend compdb log(const compdb &z) { return compdb(log(norm(z))*0.5, arg(z)); }
    friend compdb coth(const compdb &z);
    friend compdb atan(const compdb &z);
};

namespace Maths
{
  extern const compdb ci(0.0, 1.0);
  extern const compdb cz(0.0, 0.0);
}

inline double abs(const compdb &z)
{
  double R = std::fabs(z.re), I = std::fabs(z.im);
  return R > I ? R*sqrt(1.0+(I/R)*(I/R)) : ( I == 0.0 ? 0.0 : I*sqrt(1.0+(R/I)*(R/I)) );
}

inline compdb compdb::operator-() const
{
  return compdb(-re,-im);
}

inline compdb& compdb::operator= (const compdb &z)
{
  if (this == &z)
    return *this;
  re = z.re;
  im = z.im;
  return *this;
}

inline compdb& compdb::operator= (const double r)
{
  re = r;
  im = 0.0;
  return *this;
}

inline compdb& compdb::operator+= (const compdb &z)
{
  re += z.re;
  im += z.im;
  return *this;
}

inline compdb& compdb::operator+= (const double r)
{
  re += r;
  return *this;
}

inline compdb& compdb::operator-= (const compdb &z)
{
  re -= z.re;
  im -= z.im;
  return *this;
}

inline compdb& compdb::operator-= (const double r)
{
  re -= r;
  return *this;
}

inline compdb& compdb::operator*= (const compdb &z)
{
  double a = re*z.re - im*z.im;
  im = re*z.im + im*z.re;
  re = a;
  return *this;
}

inline compdb& compdb::operator*= (const double r)
{
  re *= r;
  im *= r;
  return *this;
}

inline compdb& compdb::operator/= (const compdb &z)
{
  //double norm_ = norm(z);
  //double a = re*z.re + im*z.im;
  //im = im*z.re - re*z.im;
  //re = a/norm_;
  //im /= norm_;
  //return *this;
  if (std::fabs(z.re) < std::fabs(z.im)) {
    double q = z.re/z.im;
    double d = z.re*q + z.im;
    double r = re;
    re = (r*q+im)/d;
    im = (im*q-r)/d;
  } else {
    double q = z.im/z.re;
    double d = z.re + z.im*q;
    double r = re;
    re = (r+im*q)/d;
    im = (im-r*q)/d;
  }
  return *this;
}

inline compdb& compdb::operator/= (const double r)
{
  re /= r;
  im /= r;
  return *this;
}

inline bool operator==  (const compdb &x, const compdb &y) { return (x.re==y.re) && (x.im==y.im); }
inline bool operator!=  (const compdb &x, const compdb &y) { return !(operator==(x,y)); }
inline compdb operator+ (const compdb &x, const compdb &y) { return compdb(x.re + y.re, x.im + y.im); }
inline compdb operator- (const compdb &x, const compdb &y) { return compdb(x.re - y.re, x.im - y.im); }
inline compdb operator* (const compdb &x, const compdb &y) { return compdb(x.re*y.re-x.im*y.im, x.re*y.im+x.im*y.re); }
inline compdb operator/ (const compdb &x, const compdb &y)
{
  //double norm_ = norm(y);
  //return compdb{(x.re*y.re+x.im*y.im)/norm_, (x.im*y.re-x.re*y.im)/norm_};
  if (std::fabs(y.re) < std::fabs(y.im)) {
    double q = y.re/y.im;
    double d = y.re*q + y.im;
    return compdb((x.re*q+x.im)/d, (x.im*q-x.re)/d);
  } else {
    double q = y.im/y.re;
    double d = y.re + y.im*q;
    return compdb((x.re+x.im*q)/d, (x.im-x.re*q)/d);
  }
}

inline bool operator==  (const compdb &x, double y) { return (x.re==y) && (x.im==0); }
inline bool operator==  (double x, const compdb &y) { return (x==y.re) && (y.im==0); }
inline bool operator!=  (const compdb &x, double y) { return !(operator==(x,y)); }
inline bool operator!=  (double x, const compdb &y) { return !(operator==(x,y)); }
inline compdb operator+ (const compdb &x, double y) { return compdb(x.re + y, x.im); } 
inline compdb operator+ (double x, const compdb &y) { return compdb(x + y.re, y.im); }
inline compdb operator- (const compdb &x, double y) { return compdb(x.re - y, x.im); }
inline compdb operator- (double x, const compdb &y) { return compdb(x - y.re,-y.im); }
inline compdb operator* (const compdb &x, double y) { return compdb(x.re*y, x.im*y); }
inline compdb operator* (double x, const compdb &y) { return compdb(x*y.re, x*y.im); }
inline compdb operator/ (const compdb &x, double y) { return compdb(x.re/y, x.im/y); }
inline compdb operator/ (double x, const compdb &y)
{
  //double xnorm = x/norm(y);
  //return compdb{xnorm*y.re, -xnorm*y.im};
  if (std::fabs(y.re) < std::fabs(y.im)) {
    double q = y.re/y.im;
    double d = y.re*q + y.im;
    return compdb(x*q/d, -x/d);
  } else {
    double q = y.im/y.re;
    double d = y.re + y.im*q;
    return compdb(x/d, -x*q/d);
  }
}

inline compdb polar(double r, double phi)
{
  if (r == 0.0 && !std::isnan(phi))
    //return compdb{0.0,0.0};
    return 0.0;
  else
    return compdb(r*cos(phi), r*sin(phi));
}


inline compdb sqrt(const compdb &z1, const compdb &z2)
{
  double phi1 = atan2(z1.im, z1.re);
  double phi2 = atan2(z2.im, z2.re);
  double phi = 0.5*(phi1+phi2);
  double r = sqrt(sqrt( norm(z1)*norm(z2) ));
  return compdb(r*cos(phi),r*sin(phi));
}

inline compdb sqrt(const compdb &x)
{
  double r = abs(x);
  double nr, ni;
  if (r == 0.0) {
    nr = ni = r;
  } else if (real(x) > 0) {
    nr = sqrt(0.5*(r + real(x)));
    ni = imag(x)/nr/2;
  } else {
    ni = sqrt(0.5*(r - real(x)));
    if (imag(x) < 0)
      ni = -ni;
    nr = imag(x)/ni/2;
  }
  return compdb(nr, ni);
}

inline compdb sqrt_(const compdb &x)
{
  double r = sqrt(abs(x));
  double phi = atan2(x.im, x.re);
  return compdb(r*cos(0.5*phi), r*sin(0.5*phi));
}

inline compdb exp(const compdb &z)
{
  double er = exp(z.re);
  if (z.im == 0.0) return compdb(er,0.0);
  if (er == 0.0) return compdb(0.0,0.0);

  if (std::isinf(er)) {
    if (std::isnan(z.im) || std::isinf(z.im))
      return compdb(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN());
    else
      return compdb(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
  }
  return compdb(er*cos(z.im), er*sin(z.im));
}

inline compdb coth(const compdb &z)
{
  double ex = exp(z.re), emx = exp(-z.re);
  double cy = cos(z.im), sy = sin(z.im);
  double u1 = ex+emx, u2 = ex-emx;
  double i = 1./( (u2*cy)*(u2*cy)+(u1*sy)*(u1*sy) );
  return compdb(u1*u2*i, cy*sy*(u2*u2-u1*u1)*i);
}

inline compdb atan(const compdb &z)
{
  double x = z.re, y = z.im, x2 = x*x, y2=y*y;
  double den = 1./(x2+(y+1.)*(y+1.));
  double r = 0.5*atan2(2*x*den, (1-x2-y2)*den);
  double i = -0.25*log(den*(x2+(y-1.)*(y-1.)));
  return compdb(r,i);
}

inline std::ostream& operator<< (std::ostream &oS, const compdb &z)
{
  int w = oS.width();
  oS << std::setw(w)<< z.re << "\t" << std::setw(w) << z.im;
  return oS;
}

inline std::istream& operator>> (std::istream &iS, compdb &z){
  iS >> z.re >> z.im;
  return iS;
}


#endif /* end of include guard: COMPLEX_H */
