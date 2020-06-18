#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <numeric>
#include "assert.h"

template <typename T> std::vector<T> operator+(const std::vector<T> &v, const std::vector<T> &w);
template <typename T> std::vector<T> operator+(const std::vector<T> &v, T a);
template <typename T> std::vector<T> operator+(T a, const std::vector<T> &v);
template <typename T> std::vector<T> operator*(const std::vector<T> &v, const std::vector<T> &w);
template <typename T> std::vector<T> operator*(const std::vector<T> &v, T a);
template <typename T> std::vector<T> operator*(T a, const std::vector<T> &v);

template <typename T> T sum(const std::vector<T> &v);
template <typename T> T product(const std::vector<T> &v);
template <typename T> T dot(const std::vector<T> &v, const std::vector<T> &w);
template <typename T> std::vector<T> cross(const std::vector<T> &v, const std::vector<T> &w);


template <typename T> std::vector<T> operator+(const std::vector<T> &v, const std::vector<T> &w)
{
  assert0(v.size() == w.size());
  std::vector<T> vec {v};
  for (int i = 0; i < v.size(); i++)
    vec[i] += w[i];
  return vec;
}

template <typename T> std::vector<T> operator+(const std::vector<T> &v, T a)
{
  std::vector<T> vec {v};
  for (auto &i : vec)
    i += a;
  return vec;
}

template <typename T> std::vector<T> operator+(T a, const std::vector<T> &v) { return v + a; }

template <typename T> std::vector<T> operator*(const std::vector<T> &v, const std::vector<T> &w)
{
  assert0(v.size() == w.size());
  std::vector<T> vec {v};
  for (int i = 0; i < v.size(); i++)
    vec[i] *= w[i];
  return vec;
}

template <typename T> std::vector<T> operator*(const std::vector<T> &v, T a)
{
  std::vector<T> vec {v};
  for (auto &i : vec) i *= a;
  return vec;
}

template <typename T> std::vector<T> operator*(T a, const std::vector<T> &v) { return v * a; }

template <typename T> T sum(const std::vector<T> &v)
{
  return std::accumulate(v.begin(), v.end(), T(0));
}

template <typename T> T product(const std::vector<T> &v)
{
  return std::accumulate(v.begin(), v.end(), T(1), std::multiplies<T>());
}

template <typename T> T dot(const std::vector<T> &v, const std::vector<T> &w)
{
  assert0(v.size() == w.size());
  return std::inner_product(v.begin(), v.end(), w.begin(), T(0));
}

template <typename T> std::vector<T> cross(const std::vector<T> &v, const std::vector<T> &w)
{
  assert( (v.size() == 3 && v.size() == w.size()), "Cross product: only 3d!");
  return std::vector<T>{v[1]*w[2] - v[2]*w[1], v[2]*w[0] - v[0]*w[2], v[0]*w[1] - v[1]*w[0]};
}

#endif /* end of include guard VECTOR_H */
