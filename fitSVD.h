#ifndef FITSVD_H
#define FITSVD_H

#include <iostream>
#include <iomanip>
#include "maths.h"
#include "function.h"
#include "mesh.h"
#include "gamma.h"
#include "svd.h"

class FitSVD
{ // Object for general linear least-squares fitting using singular value decomposition.
  // the $\chi^2$ minimization that fits for the coefficients a[0..M-1] of a function that depends linearly on a,
  // $y = \sum_i a_i f_i(x)$.
  //    Example is given at the bottom
  private:
    int N;         // number of data points
    int M;         // number of fitting coefficients

  public:
    F1d<double> a; // the fitting coefficients a[M].
    F2d<double> C; // the covariance matrix C[M][M],
    double chi2;   // $\chi^2$
    double q;      // goodness-of-fit probability Q(chi2,N-M) = 1- cdf of chi2 distribution = incomplete Gamma function
    double s;      // the estimated error of each point.

  private:
    F2d<double> A; // design matrix in "$\vec{a} s.t. min(\chi^2 = |A \cdot \vec{a} - \vec{b}|^2)".
    F1d<double> b; // vector b in "$\vec{a} s.t. min(\chi^2 = |A \cdot \vec{a} - \vec{b}|^2)".

    double tol;
    bool existS;

  public:
    template <class functor> // A functor f(double) = F1d<double> containing M basis functions evaluated at x.
    FitSVD(const M1d &x, const F1d<double> &y, const F1d<double> &sig, functor &f, const double TOL = 1e-12);
    template <class functor>
    FitSVD(const M1d &x, const F1d<double> &y, functor &f, const double TOL = 1.e-12);

    ~FitSVD() { }

    void fit();
    friend std::ostream& operator<<(std::ostream &oS, const FitSVD &f);

};

template <class functor>
FitSVD::FitSVD(const M1d &x, const F1d<double> &y, const F1d<double> &sig, functor &f, const double TOL)
  : N{x.size()}, M{f(x[0]).size()}, a{M}, C{M,M}, A{N,M}, b{N}, tol{TOL}, existS{true}
{ // If TOL > 0, it is the thresh (relative to the largest singular value) for discarding small singular values.
  // If TOL <= 0, the default value in SVD is used.
  F1d<double> basis{M};
  for (int i = 0; i < N; i++) {
    basis = f(x[i]);
    for (int j = 0; j < M; j++)
      A(i,j) = basis[j]/sig[i];
    b[i] = y[i]/sig[i];
  }
}

template <class functor>
FitSVD::FitSVD(const M1d &x, const F1d<double> &y, functor &f, const double TOL)
  : N{x.size()}, M{f(x[0]).size()}, a{M}, C{M,M}, A{N,M}, b{N}, tol{TOL}, existS{false}
{ 
  F1d<double> basis{M};
  for (int i = 0; i < N; i++) {
    basis = f(x[i]);
    for (int j = 0; j < M; j++)
      A(i,j) = basis[j];
    b[i] = y[i];
  }
}

void FitSVD::fit()
{
  SVD svd(A);
  double thresh = (tol > 0. ? tol*svd.w[0] : -1.);
  svd.solve(b, a, thresh);

  chi2 = 0.0;
  for (int i = 0; i < N; i++) {
    double sum = 0.0;
    for (int j = 0; j < M; j++)
      sum += A(i,j)*a[j];
    chi2 += Maths::sqr(sum-b[i]);
  }

  int nu = N-M;
  if (existS) {
    s = 0.0;
    q = (nu > 0) ? Gamma::Q(0.5*nu,0.5*chi2) : 1.0;
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < i+1; j++) {
        double sum = 0.0;
        for (int k = 0; k < M; k++) {
          if (svd.w[k] > svd.tsh)
            sum += svd.v(i,k)*svd.v(j,k)/Maths::sqr(svd.w[k]);
        }
        C(j,i) = C(i,j) = sum;
      }
    }
  } else {
    s = chi2/nu;
    q = 1.0;
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < i+1; j++) {
        double sum = 0.0;
        for (int k = 0; k < M; k++) {
          if (svd.w[k] > svd.tsh)
            sum += svd.v(i,k)*svd.v(j,k)/Maths::sqr(svd.w[k]);
        }
        C(j,i) = C(i,j) = s*sum;
      }
    }
  }
}

inline std::ostream& operator<< (std::ostream &oS, const FitSVD &f)
{
  int w = oS.width();

  oS << "FitSVD chi2=" << std::setw(w) << f.chi2 << " existS="  << (f.existS ? "T" : "F");
  oS << " q = " << std::setw(w) << f.q << ", s = " << std::setw(w) << f.s << " (tol=" << std::setw(w) << f.tol << ")";

  return oS;
}

#ifdef EXAMPLE_123123123
class Polynomial
{ // ploynomial of degree n-1 functor as a basis function 
  private:
    int n;
  public:
    Polynomial(int n_) : n{n_+1} { }
    ~Polynomial() { }
    F1d<double> operator()(const double x)
    { 
      F1d<double> p{n};
      p[0] = 1.0;
      for (int i = 1; i < n; i++) {
        p[i] = p[i-1]*x;
      }
      return p; // = {1.0, x, x^2, ..., x^n-1}
    }
};

int main(void)
{
  ...

  // In a code, M1d x, F1d<double> data are given.
  Polynomial functor{M};
  double tol = 1.0e-6;
  FitSVD fitted(x, data, functor, tol);
  fitted.fit();

  ...
}
#endif


#endif /* end of include guard: FITSVD_H */
