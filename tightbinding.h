#ifndef TIGHTBINDING_H
#define TIGHTBINDING_H

#include <map>
#include "maths.h"
#include "complex.h"
#include "vector.h"
#include "matrix.h"
#include "symmetry.h"

class TightBinding
{ // Given hopping t matrix (hopT), unit vector A, kmesh for each axis, and positions of orbitals in terms of A,
  // calculates H_k by Fourier transformation.
  public:
    std::vector<Mat<compdb>> Hk;
    std::vector<double> weight;

  private:
    int dim;
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> B;
    std::vector<std::vector<double>> K;

  public:
    TightBinding() { }
    TightBinding(const std::vector<std::vector<double>> &A, int spaceGroup, const std::vector<int> &kmesh) { init(A, spaceGroup, kmesh); }
    TightBinding(const TightBinding &tb) : Hk(tb.Hk), dim(tb.dim), A(tb.A), B(tb.B), K(tb.K) { }
    ~TightBinding() { }

    TightBinding& operator= (const TightBinding &tb);
    int size() const { return Hk.size(); }
    void init(const std::vector<std::vector<double>> &A_, int spaceGroup, const std::vector<int> &kmesh);
    void init(const std::vector<std::vector<double>> &A_, const std::vector<std::vector<double>> &kpaths, int nkpath);
    void computeHk(const std::map<std::vector<int>,std::vector<std::vector<double>>> &hopT, const std::vector<std::vector<double>> &orbitals);

  private:
    void genReciprocal();
    void genKpoints(Symmetry &sym, const std::vector<int> &kmesh);
    void genKpaths(const std::vector<std::vector<double>> &kpaths, int nkpath);
};

TightBinding& TightBinding::operator= (const TightBinding &tb)
{
  if (this == &tb)
    return *this;
  this->Hk = tb.Hk;
  this->dim = tb.dim;
  this->A = tb.A;
  this->B = tb.B;
  this->K = tb.K;

  return *this;
}

void TightBinding::init(const std::vector<std::vector<double>> &A_, int spaceGroup, const std::vector<int> &kmesh)
{
  std::clog << "Constructing TB.\n";
  dim = A_.size();
  A = A_;

  B.resize(dim);
  for (int d = 0; d < dim; d++)
    B[d].resize(dim);
  genReciprocal();
  Symmetry sym(dim, spaceGroup);
  genKpoints(sym, kmesh);
}

void TightBinding::init(const std::vector<std::vector<double>> &A_, const std::vector<std::vector<double>> &kpaths, int nkpath)
{
  dim = A_.size();
  A = A_;

  B.resize(dim);
  for (int d = 0; d < dim; d++)
    B[d].resize(dim);
  genReciprocal();
  genKpaths(kpaths, nkpath);
}

void TightBinding::genReciprocal()
{
  for (int d = dim; d < 3; d++) {
    for (auto &a : A)
      a.push_back(0.);
    std::vector<double> a(d+1, 0.); a[d] = 1.;
    A.push_back(a);
  }
  double fac = 2.0*Maths::pi/dot(cross(A[0],A[1]),A[2]);
  for (int i = 0; i < dim; i++) {
    B[i] = fac*cross(A[(i+1)%3], A[(i+2)%3]);
    B[i].resize(dim);
  }
  for (int d = 3; d --> dim;) {
    A.pop_back();
    for (auto &a : A)
      a.pop_back();
  }
  std::clog << "Reciprocal lattice:\n";
  for (auto &bi : B) {
    std::clog << "   ";
    for (auto &bij : bi)
      std::clog << std::setw(15) << bij << ' ';
    std::clog << "\n";
  }
}

void TightBinding::genKpoints(Symmetry &sym, const std::vector<int> &kmesh)
{
  std::vector<int> n(dim);
  auto index2point = [&](size_t I, std::vector<double> &p)
  { // \vec{p} = B \vec{m}, m_i = n_i/N_i - 1/2, I = (n_1*N_2 + n_2)*N_3 + n_3 for d=3.
    n[dim-1] = I % kmesh[dim-1];
    for (int i = dim-1; i-->0;)
      n[i] = (I/=kmesh[i+1])%kmesh[i];
    for (int i = 0; i < dim; i++) {
      p[i] = 0.0;
      for (int j = 0; j < dim; j++)
        p[i] += B[j][i]*(static_cast<double>(n[j])/kmesh[j] - 0.5);
    }
  };
  auto point2index = [&](const std::vector<double> &p, size_t &I)
  { // \vec{m} = \frac{1}{2\pi} A^T \vec{p}
    for (int i = 0; i < dim; i++)
      if ( (n[i] = std::round(0.5*kmesh[i]*(dot(A[i],p)/Maths::pi + 1.0))) >= kmesh[i] )
        return false;
    I = n[0];
    for (int i = 1; i < dim; i++)
      I = kmesh[i]*I + n[i];
    return true;
  };

  std::vector<bool> fullK(product(kmesh), true); // faster than vector<int> in this case
  double fullksize = static_cast<double>(fullK.size());
  size_t estimateK = std::ceil(1.2*fullksize/sym.order());
  K.reserve(estimateK);
  weight.reserve(estimateK);

  std::vector<double> p(dim);
  for (size_t I = 0, J; I < fullK.size(); I++) {
    if (fullK[I]) {
      index2point(I,p);
      sym.findSamePoints(p);
      int w = 0;
      for (auto &sp : sym.samePoints) {
        if (point2index(sp,J) && fullK[J]) {
          fullK[J] = false;
          w++;
        }
      }
      fullK[I] = true;
      K.push_back(p);
      weight.push_back(w/fullksize);
    }
  }
  std::clog << "k-points generated: " << fullK.size() << " -> " << K.size() << " (1/" << fullksize/K.size() << ")" << std::endl;
}

void TightBinding::genKpaths(const std::vector<std::vector<double>> &kpaths, int nkpath)
{
  std::vector<std::vector<double>> highsym(kpaths.size());
  for (int l = 0; l < kpaths.size(); l++) {
    highsym[l].resize(dim);
    for (int i = 0; i < dim; i++) {
      highsym[l][i] = 0.0;
      for (int j = 0; j < dim; j++)
        highsym[l][i] += kpaths[l][j]*B[j][i];
    }
  }
  std::vector<double> len(kpaths.size()-1);
  double totlen = 0.0;
  for (int l = 0; l < kpaths.size()-1; l++) {
    len[l] = 0.0;
    for (int i = 0; i < dim; i++)
      len[l] += Maths::sqr(highsym[l+1][i] - highsym[l][i]);
    len[l] = std::sqrt(len[l]);
    totlen += len[l];
  }
  for (auto &n : len)
    n = std::round(nkpath*n/totlen);

  K.reserve(std::ceil(1.1*nkpath));
  std::vector<double> p(dim), d(dim);
  for (int l = 0; l < kpaths.size()-1; l++) {
    for (int j = 0; j < dim; j++)
      d[j] = (highsym[l+1][j] - highsym[l][j])/len[l];
    for (int i = 0; i < static_cast<int>(len[l]); i++) {
      for (int j = 0; j < dim; j++)
        p[j] = highsym[l][j]+i*d[j];
      K.push_back(p);
    }
  }
  K.push_back(highsym.back());
  std::clog << "k-points generated: " << K.size() << '\n';

  int N = 0;
  std::clog << "kpaths indices: " << N;
  for (auto &n : len)
    std::clog << ", " << (N+=static_cast<int>(n));
  std::clog << std::endl;
}

void TightBinding::computeHk(const std::map<std::vector<int>,std::vector<std::vector<double>>> &hopT, const std::vector<std::vector<double>> &orbitals)
{
  int Nb = orbitals.size();
  Hk.resize(K.size());
  for (auto &hk : Hk)
    hk.resize(Nb);

  std::vector<double> delta(dim);
  for (int k = 0; k < K.size(); k++) {
    for (int b = 0; b < Nb; b++) {
      for (int c = 0; c < Nb; c++) {
        compdb ekbc = Maths::cz;
        Hk[k](b,c) = Maths::cz;
        for (auto &t : hopT) {
          for (int i = 0; i < dim; i++)
            delta[i] = t.first[i]+orbitals[b][i]-orbitals[c][i];

          double kr = 0.0;
          for (int i = 0; i < dim; i++) {
            double r = 0.0; // \vec{r}_i = \sum_j^D \vec{\d}_j A_{ji}
            for (int j = 0; j < dim; j++)
              r += delta[j]*A[j][i];
            kr += K[k][i]*r;
          }
          ekbc += t.second[b][c]*compdb(std::cos(kr),-std::sin(kr));

          //double kr_p = 0.0, kr_m = 0.0;
          //for (int i = 0; i < dim; i++) {
          //  for (int j = 0; j < dim; j++) {
          //    kr_p += (+t.first[i]+orbitals[b][i]-orbitals[c][i])*A[i][j]*K[k][j];
          //    kr_m += (-t.first[i]+orbitals[b][i]-orbitals[c][i])*A[i][j]*K[k][j];
          //  }
          //}
          //Hk[k](b,c) += t.second[b][c]*compdb{std::cos(kr_p),-std::sin(kr_p)}
          //              +t.second[c][b]*compdb{std::cos(kr_m),-std::sin(kr_m)};
        }
        Hk[k](b,c) = ekbc;
      }
    }
  }
  std::clog << "Hk matrix computed.\n";
}

#endif /* end of include guard: TIGHTBINDING_H */
