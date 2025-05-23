#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <iomanip>
#include <initializer_list>
#include <cmath>
#include "assert.h"
#include "lapackWrap.h"

template <typename T>
class Mat 
{ // N by N square matrix.
  private:
    T* m;
    int N;
    int N2;
  public:
    Mat(int N_ = 0) : N(N_), N2(N_*N_) { m = (N) ? new T[N2] : nullptr; }
    Mat(const Mat<T> &A) { N = A.N; N2 = A.N2; m = new T[N2]; std::copy(A.m, A.m+A.N2, m); }
    Mat(std::initializer_list<T> l);
    ~Mat() { delete[] m; m = nullptr; }

    Mat<T>& operator= (const Mat<T> &A);
    T& operator()(const int i, const int j)             { assert(0<=i && i<N && 0<=j && j<N,"i,j="<<i<<","<<j<<" N="<<N); return m[i*N+j]; }
    const T& operator()(const int i, const int j) const { assert(0<=i && i<N && 0<=j && j<N,"i,j="<<i<<","<<j<<" N="<<N); return m[i*N+j]; }
    T& operator[] (const int i)             { assert(0<=i && i<N2, "i="<<i<<" N2="<<N2); return m[i]; }
    const T& operator[] (const int i) const { assert(0<=i && i<N2, "i="<<i<<" N2="<<N2); return m[i]; }

    template <typename U> friend std::ostream& operator<<(std::ostream &oS, const Mat<U> &A);
    template <typename U> friend std::istream& operator>>(std::istream &iS, Mat<U> &A);

    int row() const { return N; }
    int col() const { return N; }
    int size() const { return N2; }
    void resize(int N_ = 0);

    T* begin() const { return m; }
    T* end() const { return m+N2; }

    void inverse(Mat<T> &Ainv) const;
    void inverse();
    void eigen(Mat<T> &QT, double* E) const;
    template <typename compdb>
    void geigen(Mat<T> &B, Mat<T> &VT, compdb *E);
};

template <typename T>
Mat<T>::Mat(std::initializer_list<T> l)
{
  N2 = l.size();
  N = static_cast<int>(std::sqrt(N2));
  m = new T[N2];
  std::copy(l.begin(), l.end(), m);
}

template <typename T>
inline Mat<T>& Mat<T>::operator= (const Mat<T> &A)
{
  if (this == &A)
    return *this;
  resize(A.N);
  std::copy(A.m, A.m+A.N2, m);
  return *this;
}

template <typename T>
inline std::ostream& operator<< (std::ostream &oS, const Mat<T> &A)
{
  int w = oS.width();
  for (int i = 0; i < A.N2; i++) 
    oS << std::setw(w) << A.m[i] << '\t';
  return oS;
}

template <typename T>
inline std::istream& operator>> (std::istream &iS, Mat<T> &A)
{
  for (int i = 0; i < A.N2; i++) 
    iS >> A.m[i];
  return iS;
}

template <typename T>
inline void Mat<T>::resize(int N_)
{
  assert0(N_ >= 0);
  if (N_ > N) {
    if (m) delete[] m;
    m = new T[N_*N_];
  }
  N = N_;
  N2 = N_*N_;
}

template <typename T>
inline void Mat<T>::inverse(Mat<T> &Ainv) const
{
  switch (N) {
    case 0: {
      break;
    }
    case 1: {
      Ainv.m[0] = 1./m[0];
      break;
    }
    case 2: {
      T detinv = 1./(m[0]*m[3] - m[1]*m[2]);
      Ainv.m[0] =  detinv*m[3];
      Ainv.m[1] = -detinv*m[1];
      Ainv.m[2] = -detinv*m[2];
      Ainv.m[3] =  detinv*m[0];
      break;
    }
    case 3: {
      T detinv = 1./(m[0]*m[4]*m[8] - m[0]*m[5]*m[7] - m[1]*m[3]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7] - m[2]*m[4]*m[6]);
      Ainv.m[0] =  detinv*(m[4]*m[8] - m[5]*m[7]);
      Ainv.m[1] = -detinv*(m[1]*m[8] - m[2]*m[7]);
      Ainv.m[2] =  detinv*(m[1]*m[5] - m[2]*m[4]);
      Ainv.m[3] = -detinv*(m[3]*m[8] - m[5]*m[6]);
      Ainv.m[4] =  detinv*(m[0]*m[8] - m[2]*m[6]);
      Ainv.m[5] = -detinv*(m[0]*m[5] - m[2]*m[3]);
      Ainv.m[6] =  detinv*(m[3]*m[7] - m[4]*m[6]);
      Ainv.m[7] = -detinv*(m[0]*m[7] - m[1]*m[6]);
      Ainv.m[8] =  detinv*(m[0]*m[4] - m[1]*m[3]);
      break;
    }
    case 4: {
      T detinv = 1/(m[ 0]*(m[ 5]*(m[10]*m[15] - m[11]*m[14]) + m[ 6]*(m[11]*m[13] - m[ 9]*m[15]) + m[ 7]*(m[ 9]*m[14] - m[10]*m[13]))
                  - m[ 1]*(m[ 4]*(m[10]*m[15] - m[11]*m[14]) + m[ 6]*(m[11]*m[12] - m[ 8]*m[15]) + m[ 7]*(m[ 8]*m[14] - m[10]*m[12]))
                  + m[ 2]*(m[ 4]*(m[ 9]*m[15] - m[11]*m[13]) + m[ 5]*(m[11]*m[12] - m[ 8]*m[15]) + m[ 7]*(m[ 8]*m[13] - m[ 9]*m[12]))
                  - m[ 3]*(m[ 4]*(m[ 9]*m[14] - m[10]*m[13]) + m[ 5]*(m[10]*m[12] - m[ 8]*m[14]) + m[ 6]*(m[ 8]*m[13] - m[ 9]*m[12])));
      Ainv.m[ 0] = detinv*(m[ 5]*(m[10]*m[15] - m[11]*m[14]) + m[ 6]*(m[11]*m[13] - m[ 9]*m[15]) + m[ 7]*(m[ 9]*m[14] - m[10]*m[13]));
      Ainv.m[ 1] = detinv*(m[ 1]*(m[11]*m[14] - m[10]*m[15]) + m[ 2]*(m[ 9]*m[15] - m[11]*m[13]) + m[ 3]*(m[10]*m[13] - m[ 9]*m[14]));
      Ainv.m[ 2] = detinv*(m[ 1]*(m[ 6]*m[15] - m[ 7]*m[14]) + m[ 2]*(m[ 7]*m[13] - m[ 5]*m[15]) + m[ 3]*(m[ 5]*m[14] - m[ 6]*m[13]));
      Ainv.m[ 3] = detinv*(m[ 1]*(m[ 7]*m[10] - m[ 6]*m[11]) + m[ 2]*(m[ 5]*m[11] - m[ 7]*m[ 9]) + m[ 3]*(m[ 6]*m[ 9] - m[ 5]*m[10]));
      Ainv.m[ 4] = detinv*(m[ 4]*(m[11]*m[14] - m[10]*m[15]) + m[ 6]*(m[ 8]*m[15] - m[11]*m[12]) + m[ 7]*(m[10]*m[12] - m[ 8]*m[14]));
      Ainv.m[ 5] = detinv*(m[ 0]*(m[10]*m[15] - m[11]*m[14]) + m[ 2]*(m[11]*m[12] - m[ 8]*m[15]) + m[ 3]*(m[ 8]*m[14] - m[10]*m[12]));
      Ainv.m[ 6] = detinv*(m[ 0]*(m[ 7]*m[14] - m[ 6]*m[15]) + m[ 2]*(m[ 4]*m[15] - m[ 7]*m[12]) + m[ 3]*(m[ 6]*m[12] - m[ 4]*m[14]));
      Ainv.m[ 7] = detinv*(m[ 0]*(m[ 6]*m[11] - m[ 7]*m[10]) + m[ 2]*(m[ 7]*m[ 8] - m[ 4]*m[11]) + m[ 3]*(m[ 4]*m[10] - m[ 6]*m[ 8]));
      Ainv.m[ 8] = detinv*(m[ 4]*(m[ 9]*m[15] - m[11]*m[13]) + m[ 5]*(m[11]*m[12] - m[ 8]*m[15]) + m[ 7]*(m[ 8]*m[13] - m[ 9]*m[12]));
      Ainv.m[ 9] = detinv*(m[ 0]*(m[11]*m[13] - m[ 9]*m[15]) + m[ 1]*(m[ 8]*m[15] - m[11]*m[12]) + m[ 3]*(m[ 9]*m[12] - m[ 8]*m[13]));
      Ainv.m[10] = detinv*(m[ 0]*(m[ 5]*m[15] - m[ 7]*m[13]) + m[ 1]*(m[ 7]*m[12] - m[ 4]*m[15]) + m[ 3]*(m[ 4]*m[13] - m[ 5]*m[12]));
      Ainv.m[11] = detinv*(m[ 0]*(m[ 7]*m[ 9] - m[ 5]*m[11]) + m[ 1]*(m[ 4]*m[11] - m[ 7]*m[ 8]) + m[ 3]*(m[ 5]*m[ 8] - m[ 4]*m[ 9]));
      Ainv.m[12] = detinv*(m[ 4]*(m[10]*m[13] - m[ 9]*m[14]) + m[ 5]*(m[ 8]*m[14] - m[10]*m[12]) + m[ 6]*(m[ 9]*m[12] - m[ 8]*m[13]));
      Ainv.m[13] = detinv*(m[ 0]*(m[ 9]*m[14] - m[10]*m[13]) + m[ 1]*(m[10]*m[12] - m[ 8]*m[14]) + m[ 2]*(m[ 8]*m[13] - m[ 9]*m[12]));
      Ainv.m[14] = detinv*(m[ 0]*(m[ 6]*m[13] - m[ 5]*m[14]) + m[ 1]*(m[ 4]*m[14] - m[ 6]*m[12]) + m[ 2]*(m[ 5]*m[12] - m[ 4]*m[13]));
      Ainv.m[15] = detinv*(m[ 0]*(m[ 5]*m[10] - m[ 6]*m[ 9]) + m[ 1]*(m[ 6]*m[ 8] - m[ 4]*m[10]) + m[ 2]*(m[ 4]*m[ 9] - m[ 5]*m[ 8]));
      break;
    }
    default: {
      Ainv = *this;
      int *IPIV = new int[Ainv.N];
      int LWORK = Ainv.N2;
      T *WORK = new T[LWORK];
      int INFO;

      LAPACK::xgetrf_(&Ainv.N,&Ainv.N,Ainv.m,&Ainv.N,IPIV,&INFO);
      if (INFO) std::cerr << "xgetrf Error: INFO = " << INFO << std::endl;
      LAPACK::xgetri_(&Ainv.N,Ainv.m,&Ainv.N,IPIV,WORK,&LWORK,&INFO);
      if (INFO) std::cerr << "xgetri Error: INFO = " << INFO << std::endl;

      delete[] IPIV;
      delete[] WORK;
      break;
    }
  }
}

template <typename T>
inline void Mat<T>::inverse()
{
  switch (N) { 
    case 0: {
      break;
    }
    case 1: {
      m[0] = 1./m[0];
      break;
    }
    case 2: {
      T detinv = 1./(m[0]*m[3] - m[1]*m[2]);
      T t = m[0]; m[0] = m[3]; m[3] = t;
      m[0] *=  detinv;
      m[1] *= -detinv;
      m[2] *= -detinv;
      m[3] *=  detinv;
      break;
    }
    case 3: {
      T detinv = 1./(m[0]*m[4]*m[8] - m[0]*m[5]*m[7] - m[1]*m[3]*m[8] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7] - m[2]*m[4]*m[6]);
      Mat<T> t(*this);
      m[0] =  detinv*(t[4]*t[8] - t[5]*t[7]);
      m[1] = -detinv*(t[1]*t[8] - t[2]*t[7]);
      m[2] =  detinv*(t[1]*t[5] - t[2]*t[4]);
      m[3] = -detinv*(t[3]*t[8] - t[5]*t[6]);
      m[4] =  detinv*(t[0]*t[8] - t[2]*t[6]);
      m[5] = -detinv*(t[0]*t[5] - t[2]*t[3]);
      m[6] =  detinv*(t[3]*t[7] - t[4]*t[6]);
      m[7] = -detinv*(t[0]*t[7] - t[1]*t[6]);
      m[8] =  detinv*(t[0]*t[4] - t[1]*t[3]);
      break;
    }
    case 4: {
      T detinv = 1/(m[ 0]*(m[ 5]*(m[10]*m[15] - m[11]*m[14]) + m[ 6]*(m[11]*m[13] - m[ 9]*m[15]) + m[ 7]*(m[ 9]*m[14] - m[10]*m[13]))
                  - m[ 1]*(m[ 4]*(m[10]*m[15] - m[11]*m[14]) + m[ 6]*(m[11]*m[12] - m[ 8]*m[15]) + m[ 7]*(m[ 8]*m[14] - m[10]*m[12]))
                  + m[ 2]*(m[ 4]*(m[ 9]*m[15] - m[11]*m[13]) + m[ 5]*(m[11]*m[12] - m[ 8]*m[15]) + m[ 7]*(m[ 8]*m[13] - m[ 9]*m[12]))
                  - m[ 3]*(m[ 4]*(m[ 9]*m[14] - m[10]*m[13]) + m[ 5]*(m[10]*m[12] - m[ 8]*m[14]) + m[ 6]*(m[ 8]*m[13] - m[ 9]*m[12])));
      Mat<T> t(*this);
      m[ 0] = detinv*(t[ 5]*(t[10]*t[15] - t[11]*t[14]) + t[ 6]*(t[11]*t[13] - t[ 9]*t[15]) + t[ 7]*(t[ 9]*t[14] - t[10]*t[13]));
      m[ 1] = detinv*(t[ 1]*(t[11]*t[14] - t[10]*t[15]) + t[ 2]*(t[ 9]*t[15] - t[11]*t[13]) + t[ 3]*(t[10]*t[13] - t[ 9]*t[14]));
      m[ 2] = detinv*(t[ 1]*(t[ 6]*t[15] - t[ 7]*t[14]) + t[ 2]*(t[ 7]*t[13] - t[ 5]*t[15]) + t[ 3]*(t[ 5]*t[14] - t[ 6]*t[13]));
      m[ 3] = detinv*(t[ 1]*(t[ 7]*t[10] - t[ 6]*t[11]) + t[ 2]*(t[ 5]*t[11] - t[ 7]*t[ 9]) + t[ 3]*(t[ 6]*t[ 9] - t[ 5]*t[10]));
      m[ 4] = detinv*(t[ 4]*(t[11]*t[14] - t[10]*t[15]) + t[ 6]*(t[ 8]*t[15] - t[11]*t[12]) + t[ 7]*(t[10]*t[12] - t[ 8]*t[14]));
      m[ 5] = detinv*(t[ 0]*(t[10]*t[15] - t[11]*t[14]) + t[ 2]*(t[11]*t[12] - t[ 8]*t[15]) + t[ 3]*(t[ 8]*t[14] - t[10]*t[12]));
      m[ 6] = detinv*(t[ 0]*(t[ 7]*t[14] - t[ 6]*t[15]) + t[ 2]*(t[ 4]*t[15] - t[ 7]*t[12]) + t[ 3]*(t[ 6]*t[12] - t[ 4]*t[14]));
      m[ 7] = detinv*(t[ 0]*(t[ 6]*t[11] - t[ 7]*t[10]) + t[ 2]*(t[ 7]*t[ 8] - t[ 4]*t[11]) + t[ 3]*(t[ 4]*t[10] - t[ 6]*t[ 8]));
      m[ 8] = detinv*(t[ 4]*(t[ 9]*t[15] - t[11]*t[13]) + t[ 5]*(t[11]*t[12] - t[ 8]*t[15]) + t[ 7]*(t[ 8]*t[13] - t[ 9]*t[12]));
      m[ 9] = detinv*(t[ 0]*(t[11]*t[13] - t[ 9]*t[15]) + t[ 1]*(t[ 8]*t[15] - t[11]*t[12]) + t[ 3]*(t[ 9]*t[12] - t[ 8]*t[13]));
      m[10] = detinv*(t[ 0]*(t[ 5]*t[15] - t[ 7]*t[13]) + t[ 1]*(t[ 7]*t[12] - t[ 4]*t[15]) + t[ 3]*(t[ 4]*t[13] - t[ 5]*t[12]));
      m[11] = detinv*(t[ 0]*(t[ 7]*t[ 9] - t[ 5]*t[11]) + t[ 1]*(t[ 4]*t[11] - t[ 7]*t[ 8]) + t[ 3]*(t[ 5]*t[ 8] - t[ 4]*t[ 9]));
      m[12] = detinv*(t[ 4]*(t[10]*t[13] - t[ 9]*t[14]) + t[ 5]*(t[ 8]*t[14] - t[10]*t[12]) + t[ 6]*(t[ 9]*t[12] - t[ 8]*t[13]));
      m[13] = detinv*(t[ 0]*(t[ 9]*t[14] - t[10]*t[13]) + t[ 1]*(t[10]*t[12] - t[ 8]*t[14]) + t[ 2]*(t[ 8]*t[13] - t[ 9]*t[12]));
      m[14] = detinv*(t[ 0]*(t[ 6]*t[13] - t[ 5]*t[14]) + t[ 1]*(t[ 4]*t[14] - t[ 6]*t[12]) + t[ 2]*(t[ 5]*t[12] - t[ 4]*t[13]));
      m[15] = detinv*(t[ 0]*(t[ 5]*t[10] - t[ 6]*t[ 9]) + t[ 1]*(t[ 6]*t[ 8] - t[ 4]*t[10]) + t[ 2]*(t[ 4]*t[ 9] - t[ 5]*t[ 8]));
      break;
    }
    default: {
      int *IPIV = new int[N];
      int LWORK = N2;
      T *WORK = new T[LWORK];
      int INFO;

      LAPACK::xgetrf_(&N,&N,m,&N,IPIV,&INFO);
      if (INFO) std::cerr << "xgetrf Error: INFO = " << INFO << std::endl;
      LAPACK::xgetri_(&N,m,&N,IPIV,WORK,&LWORK,&INFO);
      if (INFO) std::cerr << "xgetri Error: INFO = " << INFO << std::endl;

      delete[] IPIV;
      delete[] WORK;
      break;
    }
  }
}

template <>
inline void Mat<double>::eigen(Mat<double> &QT, double* E) const
{
  QT = *this;
  char JOBZ='V';
  char UPLO='L';
  int LWORK = 3*N-1;
  double *WORK = new double[LWORK];
  int INFO;
  int N_ = N;

  LAPACK::dsyev_(&JOBZ, &UPLO, &N_, QT.m, &N_, E, WORK, &LWORK, &INFO);

  delete[] WORK;
  if (INFO) std::cerr << "dsyev Error: INFO = " << INFO << std::endl;
}

template <typename T>
inline void Mat<T>::eigen(Mat<T> &QT, double* E) const
{
  QT = *this;
  char JOBZ='V';
  char UPLO='L';
  int LWORK = 2*N-1;
  T *WORK = new T[LWORK];
  double *RWORK = new double[3*N-2];
  int INFO;
  int N_ = N;

  LAPACK::zheev_(&JOBZ, &UPLO, &N_, QT.m, &N_, E, WORK, &LWORK, RWORK, &INFO);

  delete[] WORK;
  delete[] RWORK;
  if (INFO) std::cerr << "zheev Error: INFO = " << INFO << std::endl;
}

template <>
template <typename compdb>
void Mat<double>::geigen(Mat<double> &B, Mat<double> &VT, compdb *E)
{
  char JOBVL = 'N';
  char JOBVR = 'V';
  int LWORK = 8*N+16;
  int INFO;
  double VL;
  
  double *ALPHAR = new double[N];
  double *ALPHAI = new double[N];
  double *BETA = new double[N];
  double *WORK = new double[LWORK];

  LAPACK::dggev_(&JOBVL, &JOBVR, &N, m, &N, B.m, &N, ALPHAR, ALPHAI, BETA, &VL, &N, VT.m, &N, WORK, &LWORK, &INFO);
  if (INFO) std::cerr << "dggev Error: INFO = " << INFO << std::endl;

  for (int i = 0; i < N; i++)
    E[i] = compdb(ALPHAR[i],ALPHAI[i])/BETA[i];

  delete[] ALPHAR;
  delete[] ALPHAI;
  delete[] BETA;
  delete[] WORK;
}


#endif /* end of include guard: MATRIX_H */
