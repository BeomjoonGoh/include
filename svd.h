#ifndef SVD_H
#define SVD_H

#include "assert.h"
#include "util.h"
#include "maths.h"
#include "function.h"

class SVD
{ // Object for singular value decomposition of a matrix A, and related functions.
  // Similar to the LAPACK routine GESVD.
  private:
    const int m, n;
    const int MaxIt;
  public:
    F2d<double> u;
    F2d<double> v; 
    F1d<double> w;
    double tsh;

  public:
    SVD(const F2d<double> &A);
    virtual ~SVD() { }

    void solve(const F1d<double> &b, F1d<double> &x, double thresh = -1.);
    void solve(const F2d<double> &b, F2d<double> &x, double thresh = -1.);

    int Rank(double thresh = -1.0);
    int nullity(double thresh = -1.0);

    F2d<double> range(double thresh = -1.0);
    F2d<double> nullspace(double thresh = -1.0);

    double inv_condition() { return (w[0] <= 0. || w[n-1] <= 0.) ? 0. : w[n-1]/w[0]; }

  private:
    void decompose();
    void reorder();
    double pythag(const double a, const double b);
};

SVD::SVD(const F2d<double> &A) : m{A.sizeD()}, n{A.size()}, MaxIt{30}, u{A}, v{n,n}, w{n} 
{ // Constructor. The SVD computation is done by decompose, and the results are sorted by reorder.
  decompose();
  reorder();
  tsh = 0.5*std::sqrt(m+n+1.)*w[0]*Maths::absEpsilon;
}

int SVD::Rank(double thresh)
{ // Return the rank of A, after zeroing any singular values smaller than thresh.
  // If thresh is negative, a default value based on estimated roundoff is used.
  tsh = (thresh >= 0.0) ? thresh : 0.5*std::sqrt(m+n+1.)*w[0]*Maths::absEpsilon;
  int nr = 0;
  for (int j = 0; j < n; j++)
    if (w[j] > tsh) nr++;
  return nr;
}

int SVD::nullity(double thresh)
{ // Return the nullity of A, after zeroing any singular values smaller than thresh. Default value as above.
  tsh = (thresh >= 0.) ? thresh : 0.5*std::sqrt(m+n+1.)*w[0]*Maths::absEpsilon;
  int nn = 0;
  for (int j = 0; j < n; j++)
    if (w[j] <= tsh) nn++;
  return nn;
}

F2d<double> SVD::range(double thresh)
{ // Give an orthonormal basis for the range of A as the columns of a returned matrix. thresh as above.
  F2d<double> R{m,Rank(thresh)};
  int nr = 0;
  for (int j = 0; j < n; j++) {
    if (w[j] > tsh) {
      for (int i = 0; i < m; i++)
        R(i,nr) = u(i,j);
      nr++;
    }
  }
  return R;
}

F2d<double> SVD::nullspace(double thresh)
{ // Give an orthonormal basis for the nullspace of A as the columns of a returned matrix. thresh as above.
  F2d<double> N{n,nullity(thresh)};
  int nn = 0;
  for (int j = 0; j < n; j++) {
    if (w[j] <= tsh) {
      for (int jj = 0; jj < n; jj++)
        N(jj,nn) = v(jj,j);
      nn++;
    }
  }
  return N;
}

void SVD::solve(const F1d<double> &b, F1d<double> &x, double thresh)
{ // Solve Ax = b for a vector x using the Moore-Penrose inverse (pseudoinverse) of A as obtained by SVD.
  // If positive, thresh is the threshold value below which singular values are considered as zero.
  // If thresh is negative, a default based on expected roundoff error is used.
  assert(b.size() == m && x.size() == n, "b.size(),m="<<b.size()<<","<<m<<", x.size(),n="<<x.size()<<","<<n);
  tsh = (thresh >= 0. ? thresh : 0.5*std::sqrt(m+n+1.)*w[0]*Maths::absEpsilon);

  F1d<double> tmp{n};
  for (int j = 0; j < n; j++) { // Calculate U^T B.
    double sum = 0.0;
    if (w[j] > tsh) {  // Nonzero result only if w_j is nonzero.
      for (int i = 0; i < m; i++) {
        sum += u(i,j)*b[i];
      }
      sum /= w[j];
    }
    tmp[j] = sum;
  }

  for (int j = 0; j < n; j++) {
    double sum = 0.0;
    for (int jj = 0; jj < n; jj++) {
      sum += v(j,jj)*tmp[jj];
    }
    x[j] = sum;
  }
}

void SVD::solve(const F2d<double> &b, F2d<double> &x, double thresh)
{ // Solves M sets of n equations A 􏰻 X D B using the pseudoinverse of A. The right-hand sides are input as
  // b[0..n-1][0..M-1], while x[0..n-1][0..M-1] returns the solutions. thresh as above.
  assert(b.sizeD() == n && x.sizeD() == n && b.size() == x.size(),
         "b.sizeD(),x.sizeD(),n="<<b.sizeD()<<","<<x.sizeD()<<","<<n<<", b.size(),x.size()="<<b.size()<<","<<x.size());
  int M = b.size();

  F1d<double> xx{n};
  for (int j = 0; j < M; j++) { // Copy and solve each column in turn.
    for (int i = 0; i < n; i++) { xx[i] = b(i,j); }
    solve(xx,xx,thresh);
    for (int i = 0; i < n; i++) { x(i,j) = xx[i]; }
  } 
}

void SVD::decompose()
{ // Given the matrix A stored in u[0..m-1][0..n-1], this routine computes its singular value decomposition,
  // A = U W V^T and stores the results in the matrices u and v, and the vector w.
  double anorm, f, g, h, s, scale;
  F1d<double> rv1{n};
  g = scale = anorm = 0.0;
  // Householder reduction to bidiagonal form.
  int l = 0;
  for (int i = 0; i < n; i++) {
    l = i+2;
    rv1[i] = scale*g;
    g = s = scale = 0.0;
    if (i < m) {
      for (int k = i; k < m; k++) scale += Maths::abs(u(k,i));
      if (scale != 0.0) {
        for (int k = i; k < m; k++) {
          u(k,i) /= scale;
          s += u(k,i)*u(k,i);
        }
        f = u(i,i);
        g = -Maths::sign(std::sqrt(s),f);
        h = f*g-s;
        u(i,i) = f-g;
        for (int j = l-1; j < n; j++) {
          s = 0.0;
          for (int k = i; k < m; k++) s += u(k,i)*u(k,j);
          f = s/h;
          for (int k = i; k < m; k++) u(k,j) += f*u(k,i);
        }
        for (int k = i; k < m; k++) u(k,i) *= scale;
      }
    }
    w[i] = scale*g;
    g = s = scale = 0.0;
    if (i+1 <= m && i+1 != n) {
      for (int k = l-1; k < n; k++) scale += Maths::abs(u(i,k));
      if (scale != 0.0) {
        for (int k = l-1; k < n; k++) {
          u(i,k) /= scale;
          s += u(i,k)*u(i,k);
        }
        f = u(i,l-1);
        g = -Maths::sign(std::sqrt(s),f);
        h = f*g-s;
        u(i,l-1) = f-g;
        for (int k = l-1; k < n; k++) rv1[k] = u(i,k)/h;
        for (int j = l-1; j < m; j++) {
          s = 0.0;
          for (int k = l-1; k < n; k++) s += u(j,k)*u(i,k);
          for (int k = l-1; k < n; k++) u(j,k) += s*rv1[k];
        }
        for (int k = l-1; k < n; k++) u(i,k) *= scale;
      }
    }
    anorm = Maths::max(anorm,(Maths::abs(w[i])+Maths::abs(rv1[i])));
  }

  // Accumulation of right-hand transformations.
  for (int i = n-1; i >= 0; i--) {
    if (i < n-1) {
      if (g != 0.0) {
        for (int j = l; j < n; j++)
          v(j,i) = (u(i,j)/u(i,l))/g; // Double division to avoid possible underflow.
        for (int j = l; j < n; j++) {
          s = 0.0;
          for (int k = l; k < n; k++) s += u(i,k)*v(k,j);
          for (int k = l; k < n; k++) v(k,j) += s*v(k,i);
        }
      }
      for (int j = l; j < n; j++)
        v(i,j) = v(j,i) = 0.0;
    }
    v(i,i) = 1.0;
    g = rv1[i];
    l = i;
  }

  // Accumulation of left-hand transformations.
  for (int i = Maths::min(m,n)-1; i >= 0; i--) {
    l = i+1;
    g = w[i];
    for (int j = l; j < n; j++) u(i,j) = 0.0;
    if (g != 0.0) {
      g = 1.0/g;
      for (int j = l; j < n; j++) {
        s = 0.0;
        for (int k = l; k < m; k++) s += u(k,i)*u(k,j);
        f = (s/u(i,i))*g;
        for (int k = i; k < m; k++) u(k,j) += f*u(k,i);
      }
      for (int j = i; j < m; j++) u(j,i) *= g;
    } else {
      for (int j = i; j < m; j++) u(j,i)=0.0;
    }
    ++u(i,i);
  }

  // Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations.
  bool flag;
  int nm;
  double c, x, y, z;
  for (int k = n-1; k >= 0; k--) {
    for (int it = 0; it < MaxIt; it++) {
      flag = true;
      for (l = k; l >= 0; l--) { // Test for splitting.
        nm = l-1;
        if (l == 0 || Maths::abs(rv1[l]) <= Maths::absEpsilon*anorm) {
          flag = false;
          break;
        }
        if (Maths::abs(w[nm]) <= Maths::absEpsilon*anorm) break;
      }

      if (flag) {
        c = 0.0; // Cancellation of rv1[l], if l > 0.
        s = 1.0;
        for (int i = l; i < k+1; i++) {
          f = s*rv1[i];
          rv1[i] = c*rv1[i];
          if (Maths::abs(f) <= Maths::absEpsilon*anorm) break;
          g = w[i];
          h = pythag(f,g);
          w[i] = h;
          h = 1.0/h;
          c = g*h;
          s = -f*h;
          for (int j = 0; j < m; j++) {
            y = u(j,nm);
            z = u(j,i);
            u(j,nm) = y*c+z*s;
            u(j,i) = z*c-y*s;
          }
        }
      }
      z = w[k];
      // Convergence.
      if (l == k) {
        if (z < 0.0) { // Singular value is made nonnegative.
          w[k] = -z;
          for (int j = 0; j < n; j++)
            v(j,k) = -v(j,k);
        }
        break;
      }
      quitif(it == MaxIt-1, "no convergence in 30 svdcmp iterations");
      // Shift from bottom 2-by-2 minor.
      x = w[l];
      nm = k-1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g = pythag(f,1.0);
      f = ((x-z)*(x+z)+h*((y/(f+Maths::sign(g,f)))-h))/x;
      c = s = 1.0;
      // Next QR transformation:
      for (int j = l; j <= nm; j++) {
        int i = j+1;
        g = rv1[i];
        y = w[i];
        h = s*g;
        g = c*g;
        z = pythag(f,h);
        rv1[j] = z;
        c = f/z;
        s = h/z;
        f = x*c+g*s;
        g = g*c-x*s;
        h = y*s;
        y *= c;
        for (int jj = 0; jj < n; jj++) {
          x = v(jj,j);
          z = v(jj,i);
          v(jj,j) = x*c+z*s;
          v(jj,i) = z*c-x*s;
        }
        z = pythag(f,h);
        w[j] = z;
        if (z) { // Rotation can be arbitrary if z = 0.
          z = 1.0/z;
          c = f*z;
          s = h*z;
        }
        f = c*g+s*y;
        x = c*y-s*g;
        for (int jj = 0; jj < m; jj++) {
          y = u(jj,j);
          z = u(jj,i);
          u(jj,j) = y*c+z*s;
          u(jj,i) = z*c-y*s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
}

void SVD::reorder()
{ // Given the output of decompose, this routine sorts the singular values, and corresponding columns of u and v,
  // by decreasing magnitude. The method is Shell's sort. (The work is negligible as com- pared to that already done
  // in decompose.) Also, signs of corresponding columns are flipped so as to maximize the number of positive elements.
  F1d<double> su{m};
  F1d<double> sv{n};

  int inc = 1;
  do {
    inc *= 3;
    inc++;
  } while (inc <= n);

  do { 
    inc /= 3;
    for (int i = inc; i < n; i++) {
      double sw = w[i];
      for (int k = 0; k < m; k++) su[k] = u(k,i);
      for (int k = 0; k < n; k++) sv[k] = v(k,i);

      int j = i;
      while (w[j-inc] < sw) {
        w[j] = w[j-inc];
        for (int k = 0; k < m; k++) u(k,j) = u(k,j-inc);
        for (int k = 0; k < n; k++) v(k,j) = v(k,j-inc);
        j -= inc;
        if (j < inc) break;
      }
      w[j] = sw;

      for (int k = 0; k < m; k++) u(k,j) = su[k];
      for (int k = 0; k < n; k++) v(k,j) = sv[k];
    }

  } while (inc > 1);

  // Flip signs.
  for (int k = 0; k < n; k++) {
    int sum = 0;
    for (int i = 0; i < m; i++) { if (u(i,k) < 0.) sum++; }
    for (int j = 0; j < n; j++) { if (v(j,k) < 0.) sum++; }
    if (sum > (m+n)/2) {
      for (int i = 0; i < m; i++) u(i,k) = -u(i,k);
      for (int j = 0; j < n; j++) v(j,k) = -v(j,k);
    }
  }
}

inline double SVD::pythag(const double a, const double b)
{ // Computes (a^2 + b^2)^1/2 without destructive underflow or overflow.
  double A = Maths::abs(a);
  double B = Maths::abs(b);
  if (A > B)
    return A*std::sqrt(1.0+Maths::sqr(B/A));
  if (B == 0.0)
    return 0.0;
  return B*std::sqrt(1.0+Maths::sqr(A/B));
}

#endif /* end of include guard: SVD_H */
