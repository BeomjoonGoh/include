#ifndef CNUMBER_H
#define CNUMBER_H

#include <iostream>
#include <complex>
#include <cmath>
#include "maths.h"
#include "assert.h"

template <typename T>
class Cnum
{ // Polar coordinate representation of a complex number:
  //    z = r e^{i\theta}.
  private:
    T m_mag;                                        // Magnitude of a complex number Cnum<T>
    T m_arg;                                        // Argument  of a complex number Cnum<T>
  
  public:
    Cnum(T mag_ = 0.0, T arg_ = 0.0);               // takes care of real numbers as well
    Cnum(const Cnum<T> &z);                         // copy constructor
    Cnum(const std::complex<T> &z);
  
    virtual ~Cnum() { }
  
    void Set(T mag_, T arg_);
    void Reset();
  
    T Mag() const;
    T Arg() const;
    std::complex<T> Cart() const;
  
    Cnum<T>  operator- () const;
    bool     operator! () const;                    // true if z = {0.0, 0.0}
  
    Cnum<T>& operator=  (const Cnum<T> &z);
    Cnum<T>& operator+= (const Cnum<T> &z);
    Cnum<T>& operator-= (const Cnum<T> &z);
    Cnum<T>& operator*= (const Cnum<T> &z);
    Cnum<T>& operator/= (const Cnum<T> &z);
  
  
    template <typename U> friend U real(const Cnum<U> &z);
    // Or this (commenting out definition far below):
    //friend T real(const Cnum<T> &z) { return z.Cart().real(); }
    template <typename U> friend U imag(const Cnum<U> &z);
    template <typename U> friend Cnum<U> conj(const Cnum<U> &z);
  
    friend Cnum<T> operator+ (const Cnum<T> &z1, const Cnum<T> &z2);
    friend Cnum<T> operator+ (const Cnum<T> &z, T x);                 // Actually these are not needed
    friend Cnum<T> operator+ (T x, const Cnum<T> &z);                 // since Cnum will be constructed
    friend Cnum<T> operator- (const Cnum<T> &z1, const Cnum<T> &z2);  // from T x and (Cnum + Cnum)
    friend Cnum<T> operator- (const Cnum<T> &z, T x);                 // will be called.
    friend Cnum<T> operator- (T x, const Cnum<T> &z);
    friend Cnum<T> operator* (const Cnum<T> &z1, const Cnum<T> &z2);
    friend Cnum<T> operator* (const Cnum<T> &z, T x);
    friend Cnum<T> operator* (T x, const Cnum<T> &z);
    friend Cnum<T> operator/ (const Cnum<T> &z1, const Cnum<T> &z2);
    friend Cnum<T> operator/ (const Cnum<T> &z, T x);
    friend Cnum<T> operator/ (T x, const Cnum<T> &z);
  
    // same problem as in real() above.
    template <typename U> friend std::ostream& operator<< (std::ostream &oS, const Cnum<U> &z);
    template <typename U> friend std::istream& operator>> (std::istream &iS, Cnum<U> &z);
  
    friend bool operator== (const Cnum<T> &z1, const Cnum<T> &z2);
    friend bool operator== (const Cnum<T> &z, T x);
    friend bool operator== (T x, const Cnum<T> &z);
    friend bool operator!= (const Cnum<T> &z1, const Cnum<T> &z2);
    friend bool operator!= (const Cnum<T> &z, T x);
    friend bool operator!= (T x, const Cnum<T> &z);
  
  private:
    void makePrin(T &arg_);

  public:
    static const Cnum<T> cni;
};


// =====  Constructors
template <typename T> Cnum<T>::Cnum(T mag_, T arg_) { Set(mag_, arg_); }
template <typename T> Cnum<T>::Cnum(const Cnum<T> &z)    : m_mag{z.m_mag}, m_arg{z.m_arg} { }
template <typename T> Cnum<T>::Cnum(const std::complex<T> &z) : m_mag{std::abs(z)}, m_arg{std::arg(z)} { }


// =====  Member functions
// Setting functions
template <typename T>
void Cnum<T>::Set(T mag_, T arg_)
{
  if (mag_ < 0.0) {
    mag_ *= -1;
    arg_ += Maths::pi;
  }
  makePrin(arg_);

  m_mag = mag_;
  m_arg = arg_;

}
template <typename T> void Cnum<T>::Reset() { Set(0.0, 0.0); }

// Get value functions
template <typename T> T Cnum<T>::Mag() const { return m_mag; }
template <typename T> T Cnum<T>::Arg() const { return m_arg; }
template <typename T> std::complex<T> Cnum<T>::Cart() const { return std::polar(m_mag, m_arg); }

// Operator overloadings
template <typename T>
Cnum<T> Cnum<T>::operator-() const
{
  T arg_ = m_arg + Maths::pi;
  makePrin(arg_);
  return Cnum<T>{m_mag, arg_};
}
template <typename T>
bool Cnum<T>::operator!() const
{
  return ( Maths::isEqual(m_mag, 0.0) && Maths::isEqual(m_arg, 0.0) );
}

template <typename T>
Cnum<T>& Cnum<T>::operator= (const Cnum<T> &z)
{
  if (this == &z) {
    return *this;
  }
  m_mag = z.m_mag;
  m_arg = z.m_arg;
  return *this;
}
template <typename T>
Cnum<T>& Cnum<T>::operator+= (const Cnum<T> &z)
{
  std::complex<T> cz = std::polar(m_mag, m_arg) + std::polar(z.m_mag, z.m_arg);
  m_mag = std::abs(cz);
  m_arg = std::arg(cz);
  return *this;
}
template <typename T>
Cnum<T>& Cnum<T>::operator-= (const Cnum<T> &z)
{
  std::complex<T> cz = std::polar(m_mag, m_arg) - std::polar(z.m_mag, z.m_arg);
  m_mag = std::abs(cz);
  m_arg = std::arg(cz);
  return *this;
}
template <typename T>
Cnum<T>& Cnum<T>::operator*= (const Cnum<T> &z)
{
  m_mag *= z.m_mag;
  m_arg += z.m_arg;
  makePrin(m_arg);
  return *this;
}
template <typename T>
Cnum<T>& Cnum<T>::operator/= (const Cnum<T> &z)
{
  assert( !Maths::isEqual(z.m_mag, 0.0) );
  m_mag /= z.m_mag;
  m_arg -= z.m_arg;
  makePrin(m_arg);
  return *this;
}


// =====  Friend functions
// Utility functions
template <typename T> T real(const Cnum<T> &z) { return z.Cart().real(); }
template <typename T> T imag(const Cnum<T> &z) { return z.Cart().imag(); }
template <typename T> Cnum<T> conj(const Cnum<T> &z) { return Cnum<T>{z.m_mag, -z.m_arg}; }

// Operator overloadings
//  Arithmetic operators
template <typename T>
Cnum<T> operator+ (const Cnum<T> &z1, const Cnum<T> &z2)
{
  return Cnum<T>{z1.Cart() + z2.Cart()};
}
template <typename T> Cnum<T> operator+ (const Cnum<T> &z, T x) { return Cnum<T>{z.Cart() + x}; }
template <typename T> Cnum<T> operator+ (T x, const Cnum<T> &z) { return z + x; }

template <typename T>
Cnum<T> operator- (const Cnum<T> &z1, const Cnum<T> &z2)
{
  return Cnum<T>{z1.Cart() - z2.Cart()};
}
template <typename T> Cnum<T> operator- (const Cnum<T> &z, T x) { return Cnum<T>{z.Cart() - x}; }
template <typename T> Cnum<T> operator- (T x, const Cnum<T> &z) { return z - x; }

template <typename T>
Cnum<T> operator* (const Cnum<T> &z1, const Cnum<T> &z2)
{
  T arg_ = z1.m_arg + z2.m_arg;
  makePrin(arg_);
  return Cnum<T>{z1.m_mag * z2.m_mag, arg_};
}
template <typename T> Cnum<T> operator* (const Cnum<T> &z, T x) { return Cnum<T>{z.m_mag * x, z.m_arg}; }
template <typename T> Cnum<T> operator* (T x, const Cnum<T> &z) { return z * x; }

template <typename T>
Cnum<T> operator/ (const Cnum<T> &z1, const Cnum<T> &z2)
{
  assert( !Maths::isEqual(z2.m_mag, 0.0) );
  T arg_ = z1.m_arg - z2.m_arg;
  makePrin(arg_);
  return Cnum<T>{z1.m_mag / z2.m_mag, arg_};
}
template <typename T>
Cnum<T> operator/ (const Cnum<T> &z, T x)
{
  assert( !Maths::isEqual(x, 0.0) );
  return Cnum<T>{z.m_mag / x, z.m_arg};
}
template <typename T>
Cnum<T> operator/ (T x, const Cnum<T> &z)
{
  assert( !Maths::isEqual(z.m_mag, 0.0) );
  return Cnum<T>{x / z.m_mag, -z.m_arg};
}

template <typename T>
std::ostream& operator<< (std::ostream &oS, const Cnum<T> &z)
{
  oS << '(' << z.m_mag << ',' << z.m_arg << ')';
  return oS;
}
template <typename T>
std::istream& operator>> (std::istream &iS, Cnum<T> &z)
{
  iS >> z.m_mag;
  iS >> z.m_arg;
  return iS;
}

//  Comparison operators
template <typename T>
bool operator== (const Cnum<T> &z1, const Cnum<T> &z2)
{
  std::complex<T> cz1 = std::polar(z1.m_mag, z1.m_arg);
  std::complex<T> cz2 = std::polar(z2.m_mag, z2.m_arg);
  return ( Maths::isEqual(cz1.real(), cz2.real()) && Maths::isEqual(cz1.imag(), cz2.imag()) );
}
template <typename T>
bool operator== (const Cnum<T> &z, T x)
{
  std::complex<T> cz = std::polar(z.m_mag, z.m_arg);
  if ( !Maths::isEqual(cz.imag(), 0.0) ) {
    return false;
  }
  return Maths::isEqual(cz.real(), x);
}
template <typename T> bool operator== (T x, const Cnum<T> &z) { return (z == x); }
template <typename T> bool operator!= (const Cnum<T> &z1, const Cnum<T> &z2) { return !(z1 == z2); }
template <typename T> bool operator!= (const Cnum<T> &z, T x) { return !(z == x); }
template <typename T> bool operator!= (T x, const Cnum<T> &z) { return !(z == x); }


// =====  Private functions
template <typename T>
void Cnum<T>::makePrin(T &arg_)
{
  if ( abs(arg_) <= Maths::pi ) { return; }

  T twopi = 2.0*Maths::pi;
  if (arg_ < 0.0) { twopi *= -1; }

  while (!(abs(arg_) <= Maths::pi)) {
    arg_ -= twopi;
  }
}

// ===== Static Members
Cnum<float> cnif {1.0, Maths::pi/2.0};
Cnum<double> cni {1.0, Maths::pi/2.0};
Cnum<long double> cnil {1.0, Maths::pi/2.0};

//// Explicit template instantiations
//template class Cnum<float>;
//template class Cnum<double>;
//template class Cnum<long double>;

#endif /* end of include guard: CNUMBER_H */
