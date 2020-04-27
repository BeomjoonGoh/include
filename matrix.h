#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <iomanip>
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
    virtual ~Mat() { delete[] m; m = nullptr; }

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

    void inverse(Mat<T> &Ainv);
    void eigen(Mat<T> &QT, T* E) const;
};

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
inline std::istream& operator>> (std::istream &iS, Mat<T> &A){
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
inline void Mat<T>::inverse(Mat<T> &Ainv)
{
  Ainv = *this;
  int *IPIV = new int[N];
  int LWORK = N2;
  T *WORK = new T[LWORK];
  int INFO;

  LAPACK::xgetrf_(&N,&N,Ainv.m,&N,IPIV,&INFO);
  if (INFO) std::cerr << "xgetrf Error: INFO = " << INFO << std::endl;
  LAPACK::xgetri_(&N,Ainv.m,&N,IPIV,WORK,&LWORK,&INFO);
  if (INFO) std::cerr << "xgetri Error: INFO = " << INFO << std::endl;

  delete[] IPIV;
  delete[] WORK;
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
  assert(INFO == 0, "Error: dgeev returned error code " << INFO);
}

#endif /* end of include guard: MATRIX_H */
