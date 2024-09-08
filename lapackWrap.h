#ifndef LAPACKWRAP_H
#define LAPACKWRAP_H

#include "complex.h"

namespace LAPACK {
  extern "C" {
    void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
    void zgetrf_(int *M, int *N, compdb *A, int *LDA, int *IPIV, int *INFO);
    void dgetri_(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);
    void zgetri_(int *N, compdb *A, int *LDA, int *IPIV, compdb *WORK, int *LWORK, int *INFO);
    void dsyev_ (char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *info);
    void zheev_ (char *JOBZ, char *UPLO, int *N, compdb *A, int *LDA, double *W, compdb *WORK, int *LWORK, double *RWORK, int *info);
    void dggev_ (char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *B, int *LDB, double *ALPHAR, double *ALPHAI, double *BETA, double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO);
  }
  inline void xgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO) { dgetrf_(M, N, A, LDA, IPIV, INFO); }
  inline void xgetrf_(int *M, int *N, compdb *A, int *LDA, int *IPIV, int *INFO) { zgetrf_(M, N, A, LDA, IPIV, INFO); }
  inline void xgetri_(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO) { dgetri_(N, A, LDA, IPIV, WORK, LWORK, INFO); }
  inline void xgetri_(int *N, compdb *A, int *LDA, int *IPIV, compdb *WORK, int *LWORK, int *INFO) { zgetri_(N, A, LDA, IPIV, WORK, LWORK, INFO); }
} // LAPACK

#endif /* end of include guard: LAPACKWRAP_H */
