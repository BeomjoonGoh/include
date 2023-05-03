#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include "maths.h"
#include "matrix.h"

class Symmetry
{
  public:
    std::vector<std::vector<double>> samePoints;
  private:
    static const std::map<std::string, std::vector<Mat<double>>> cpg1;
    static const std::map<std::string, std::vector<Mat<double>>> cpg2;
    static const std::map<std::string, std::vector<Mat<double>>> cpg3;
    const std::vector<Mat<double>> *cpg;
    bool shrunken;

  public:
    Symmetry(int dim = 3, int spaceGroup = 1);
    ~Symmetry() { }
    void findSamePoints(const std::vector<double> &p);
    void makeSamePointsUnique();
    int order() const { return cpg->size(); }

  private:
    std::string getPointGroup1d(int spaceGroup);
    std::string getPointGroup2d(int spaceGroup);
    std::string getPointGroup3d(int spaceGroup);
};

#include "symmetry_cpg.h"

Symmetry::Symmetry(int dim, int spaceGroup)
{
  std::string pointGroup;
  switch (dim) {
    case 1:
      pointGroup = getPointGroup1d(spaceGroup);
      cpg = &cpg1.at(pointGroup);
      break;
    case 2:
      pointGroup = getPointGroup2d(spaceGroup);
      cpg = &cpg2.at(pointGroup);
      break;
    case 3:
      pointGroup = getPointGroup3d(spaceGroup);
      cpg = &cpg3.at(pointGroup);
      break;
    default:
      quitif(true, "Weird dimension, d=" << dim);
  }
  samePoints.resize(cpg->size());
  for (auto &sp : samePoints)
    sp.resize(dim);
  std::clog << "Symmetry: dim=" << dim
            << ", order of crystallographic point group '" << pointGroup << "'=" << cpg->size() << std::endl;
  shrunken = false;
}

std::string Symmetry::getPointGroup1d(int s)
{
  quitif(s <= 0 || s > 2, "Can't recognize symmetry, s=" << s);
  if (s <=   1) return "C1";
  if (s <=   2) return "D1";
  return "None";
}
std::string Symmetry::getPointGroup2d(int s)
{
  quitif(s <= 0 || s > 18, "Can't recognize symmetry, s=" << s);
  if (s <=   1) return "C1";
  if (s <=   2) return "C2";
  if (s <=   5) return "D1";
  if (s <=   9) return "D2";
  if (s <=  10) return "C4";
  if (s <=  13) return "D4";
  if (s <=  14) return "C3";
  if (s <=  16) return "D3";
  if (s <=  17) return "C6";
  if (s <=  18) return "D6";
  return "None";
}
std::string Symmetry::getPointGroup3d(int s)
{
  quitif(s <= 0 || s > 230, "Can't recognize symmetry, s=" << s);
  if (s <=   1) return "C1" ;
  if (s <=   2) return "Ci" ;
  if (s <=   5) return "C2" ;
  if (s <=   9) return "C1h";
  if (s <=  15) return "C2h";
  if (s <=  24) return "D2" ;
  if (s <=  46) return "C2v";
  if (s <=  74) return "D2h";
  if (s <=  80) return "C4" ;
  if (s <=  82) return "S4" ;
  if (s <=  88) return "C4h";
  if (s <=  98) return "D4" ;
  if (s <= 110) return "C4v";
  if (s <= 122) return "D2d";
  if (s <= 142) return "D4h";
  if (s <= 146) return "C3" ;
  if (s <= 148) return "S6" ;
  if (s <= 155) return "D3" ;
  if (s <= 161) return "C3v";
  if (s <= 167) return "D3d";
  if (s <= 173) return "C6" ;
  if (s <= 174) return "C3h";
  if (s <= 176) return "C6h";
  if (s <= 182) return "D6" ;
  if (s <= 186) return "C6v";
  if (s <= 190) return "D3h";
  if (s <= 194) return "D6h";
  if (s <= 199) return "T"  ;
  if (s <= 206) return "Th" ;
  if (s <= 214) return "O"  ;
  if (s <= 220) return "Td" ;
  if (s <= 230) return "Oh" ;
  return "None";
}

inline void Symmetry::findSamePoints(const std::vector<double> &p)
{
  if (shrunken) {
    int sizeBefore = samePoints.size();
    samePoints.resize(cpg->size());
    for (int op = sizeBefore; op < cpg->size(); op++)
      samePoints[op].resize(p.size());
    shrunken = false;
  }
  for (int op = 0; op < samePoints.size(); op++) {
    for (int i = 0; i < p.size(); i++) {
      samePoints[op][i] = 0.0;
      for (int j = 0; j < p.size(); j++)
        samePoints[op][i] += (*cpg)[op](i,j)*p[j];
    }
  }
}

void Symmetry::makeSamePointsUnique()
{
  std::sort(samePoints.begin(), samePoints.end(),
      [](const std::vector<double> &A, const std::vector<double> &B) -> bool {
        for (int i = 0; i < A.size(); i++)
          if (!Maths::isEqual(A[i], B[i], 1.0e-14))
            return A[i] < B[i];
        return false;
      });
  auto it = std::unique(samePoints.begin(), samePoints.end(),
      [](const std::vector<double> &A, const std::vector<double> &B) -> bool {
        for (int i = 0; i < A.size(); i++)
          if (!Maths::isEqual(A[i], B[i], 1.0e-14))
            return false;
        return true;
      });
  samePoints.erase(it, samePoints.end());
  shrunken = cpg->size() != samePoints.size();
}

#endif /* end of include guard: SYMMETRY_H */
