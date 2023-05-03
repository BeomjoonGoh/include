#ifndef TIGHTBINDING_H
#define TIGHTBINDING_H

#include <map>
#include "util.h"
#include "maths.h"
#include "complex.h"
#include "vector.h"

class TightBinding
{ // Given hopping t matrix (hopT), unit vector A, kmesh for each axis, and positions of orbitals in terms of A,
  // calculates H_k by Fourier transformation.

  using Vec2d = std::vector<std::vector<double>>;
  using Map = std::map<std::vector<int>,std::vector<std::vector<double>>>;

  public:
    std::vector<std::vector<std::vector<compdb>>> Hk;
  private:
    int dim;
    Vec2d A;
    Vec2d B;
    Vec2d K;
    int kSize;

  public:
    TightBinding() { }
    TightBinding(const int dim, const Vec2d &A, const std::vector<int> &kmesh) { init(dim, A, kmesh); }
    TightBinding(const TightBinding &tb) : Hk(tb.Hk), dim(tb.dim), A(tb.A), B(tb.B), K(tb.K), kSize(tb.kSize) { }
    ~TightBinding() { }

    TightBinding& operator= (const TightBinding &tb);
    double operator()(const int b, const int c, const int k) { return Hk[b][c][k].real(); }
    const double operator()(const int b, const int c, const int k) const { return Hk[b][c][k].real(); }
    int size() const { return kSize; }
    void computeHk(const int Nb, const Map &hopT, const Vec2d &orbitals);
    void init(const int dim_, const Vec2d &A_, const std::vector<int> &kmesh);

  private:
    void genReciprocal();
    int  genKpoints(const std::vector<int> &kmesh);
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
  this->kSize = tb.kSize;

  return *this;
}

void TightBinding::init(const int dim_, const Vec2d &A_, const std::vector<int> &kmesh)
{
  std::clog << "Constructing TB.\n";

  dim = dim_;
  A = A_;

  B.resize(dim);
  for (int d = 0; d < dim; d++)
    B[d].resize(dim);

  genReciprocal();
  kSize = genKpoints(kmesh);
}

void TightBinding::genReciprocal()
{
  if (dim==1) {
    B[0][0] = 2*Maths::pi/A[0][0];
  } else if (dim==2) {
    A[0].push_back(0.);
    A[1].push_back(0.);
    A.push_back(std::vector<double>{0.,0.,1.});
    double volUC = dot(cross(A[0],A[1]),A[2]);
    B[0] = cross(A[1],A[2]) * (2*Maths::pi/volUC);
    B[1] = cross(A[2],A[0]) * (2*Maths::pi/volUC);
    A.pop_back();
    for (int d = 0; d < dim; d++) {
      A[d].pop_back();
      B[d].pop_back();
    }
  } else if (dim==3) {
    double volUC = dot(cross(A[0],A[1]),A[2]);
    B[0] = cross(A[1],A[2]) * (2*Maths::pi/volUC);
    B[1] = cross(A[2],A[0]) * (2*Maths::pi/volUC);
    B[2] = cross(A[0],A[1]) * (2*Maths::pi/volUC);
  }
  std::clog << "Reciprocal lattice:\n";
  for (auto &b : B)
    std::clog << "   "  << Util::outV(b,15) << "\n";
}

int TightBinding::genKpoints(const std::vector<int> &kmesh)
{ // recursion?
  std::vector<double> fac(dim);
  if (dim==1) {
    for (int k0 = 0; k0 < kmesh[0]; k0++) {
      fac[0] = static_cast<double>(k0) / kmesh[0];
      K.push_back(B[0]*fac[0]);
    }
  } else if (dim==2) {
    for (int k0 = 0; k0 < kmesh[0]; k0++) {
      fac[0] = static_cast<double>(k0) / kmesh[0];
      for (int k1 = 0; k1 < kmesh[1]; k1++) {
        fac[1] = static_cast<double>(k1) / kmesh[1];
        K.push_back( B[0]*fac[0] + B[1]*fac[1] );
      }
    }
  } else if (dim==3) {
    for (int k0 = 0; k0 < kmesh[0]; k0++) {
      fac[0] = static_cast<double>(k0) / kmesh[0];
      for (int k1 = 0; k1 < kmesh[1]; k1++) {
        fac[1] = static_cast<double>(k1) / kmesh[1];
        for (int k2 = 0; k2 < kmesh[2]; k2++) {
          fac[2] = static_cast<double>(k2) / kmesh[2];
          K.push_back( B[0]*fac[0] + B[1]*fac[1] + B[2]*fac[2] );
        }
      }
    }
  }
  return K.size();
}

void TightBinding::computeHk(const int Nb, const Map &hopT, const Vec2d &orbitals)
{
  Hk.resize(Nb);
  for (int b = 0; b < Nb; b++) {
    Hk[b].resize(Nb);
    for (int c = 0; c < Nb; c++)
      Hk[b][c].resize(kSize);
  }

  std::vector<double> delta(dim);
  #pragma omp parallel for
  for (int k = 0; k < kSize; k++) {
    for (int b = 0; b < Nb; b++) {
      for (int c = 0; c < Nb; c++) {
        compdb ekbc{Maths::cz};
        //Hk[b][c][k] = 0.0;
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
          ekbc += t.second[b][c]*compdb{std::cos(kr),-std::sin(kr)};
          //double kr_p = 0.0, kr_m = 0.0;
          //for (int i = 0; i < dim; i++) {
          //  for (int j = 0; j < dim; j++) {
          //    kr_p += (+t.first[i]+orbitals[b][i]-orbitals[c][i])*A[i][j]*K[k][j];
          //    kr_m += (-t.first[i]+orbitals[b][i]-orbitals[c][i])*A[i][j]*K[k][j];
          //  }
          //}
          //Hk[b][c][k] += t.second[b][c]*compdb{std::cos(kr_p),-std::sin(kr_p)}
          //              +t.second[c][b]*compdb{std::cos(kr_m),-std::sin(kr_m)};
        }
        Hk[b][c][k] = ekbc;
      }
    }
  }
  std::clog << "Hk matrix computed.\n";
}

#endif /* end of include guard: TIGHTBINDING_H */
