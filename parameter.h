#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

class ParamBase
{ // For polymorphism
  public:
    ParamBase() { }
    virtual ~ParamBase() { }

    virtual bool isAssigned() const = 0;
    virtual void print(std::ostream &os) const = 0;
    virtual void read(std::istream &is, const std::string &s, char delim) = 0;
};

template <typename T>
class Param : public ParamBase
{
  private:
    const std::string name;
    T value;
    bool assigned;
  public:
    Param(const std::string &name_, const T &value_) : name(name_), value(value_), assigned(false) { }
    virtual ~Param() { }

    operator T() const { return value; } // implicit conversion

    void assign(const T &p) { value = p; }
    bool isAssigned() const { return assigned; }
    void print(std::ostream &os) const;
    void read(std::istream &is, const std::string &s, char delim);
};

class ParamSet
{
  private:
    std::vector<ParamBase*> PS;
    bool readbefore;
  public:
    ParamSet(const std::initializer_list<ParamBase*> &l) : PS(l), readbefore(false) { }
    ~ParamSet() { }

    void Log(std::ostream &os, const std::string &name);
    bool Read(const std::string &inputf, char delim = '\0');
};

template <typename T> std::ostream& operator<< (std::ostream &os, const std::vector<T> &v);
template <typename T> std::istream& operator>> (std::istream &is,       std::vector<T> &v);
template <typename T, typename U> std::ostream& operator<< (std::ostream &os, const std::map<T,U> &m);
template <typename T, typename U> std::istream& operator>> (std::istream &is,       std::map<T,U> &m);

// Param
template <typename T>
void Param<T>::print(std::ostream &os) const
{
  os.setf(std::ios::left,std::ios_base::adjustfield);
  os << std::setw(10) << name;
  os.unsetf(std::ios_base::adjustfield);
  os << " = " << value;
}

template <typename T>
void Param<T>::read(std::istream &is, const std::string &s, char delim)
{
  const int ssMax = 2000;
  if (s == name) {
    assigned = true;
    if (delim == '\0') is >> value;
    else               is.ignore(ssMax,delim) >> value;
  }
}
template <>
void Param<bool>::read(std::istream &is, const std::string &s, char delim)
{
  const int ssMax = 2000;
  if (s == name) {
    assigned = true;
    std::string sTF;
    if (delim == '\0') is >> sTF;
    else               is.ignore(ssMax,delim) >> sTF;
    std::transform(sTF.begin(), sTF.end(), sTF.begin(), ::tolower);
    std::istringstream(sTF) >> std::boolalpha >> value;
  }
}

// ParamSet
void ParamSet::Log(std::ostream &os, const std::string &name)
{
  os << "============================================================\n";
  os << "ParamSet::log() '"<< name << "'\n";
  os << std::boolalpha;
  for (auto &p : PS) {
    os << "  ";
    p->print(os);
    os << "\n";
  }
  os << std::noboolalpha;
  os << "============================================================\n";
}

bool ParamSet::Read(const std::string &inputf, char delim)
{
  if (readbefore) {
    std::cerr << "Don't read more than once!\n";
    return true;
  }
  std::ifstream inf(inputf);
  if (!inf.good()) {
    std::clog << "Missing input file! " << "No " << inputf << ". Proceeds only with default values.\n";
    return false;
  }
  readbefore = true;
  const int ssMax = 2000;
  while (inf) {
    std::string strInput;
    inf >> strInput;
    if (strInput[0] == '#') inf.ignore(ssMax,'\n');
    for (auto &p : PS) {
      if (p->isAssigned())
        continue;
      p->read(inf, strInput, delim);
    }
  }
  return true;
}

// operator <<,>> for vector and map
template <typename T>
std::ostream& operator<< (std::ostream &os, const std::vector<T> &v)
{
  os << "[ ";
  for (auto &i : v)
    os << i << ", ";
  os << "]";
  return os;
}

template <typename T, typename U>
std::ostream& operator<< (std::ostream &os, const std::map<T,U> &m)
{
  os << "{ ";
  for (auto &i : m)
    os << i.first << " : " << i.second << ", ";
  os << "}";
  return os;
}

template <typename T>
std::istream& operator>> (std::istream &is, std::vector<T> &v)
{
  char paren, comma;
  is >> std::ws >> paren;
  T val;
  v.resize(0);
  while ((is>>std::ws).peek() != ')' && is.peek() != ']') {
    is >> val;
    if ((is>>std::ws).peek() == ',') is >> comma;
    v.push_back(val);
  }
  is >> std::ws >> paren;
  return is;
}

template <typename T, typename U>
std::istream& operator>> (std::istream &is, std::map<T,U> &m)
{
  char paren, colon, comma;
  is >> std::ws >> paren;
  T key;
  U val;
  m.clear();
  while ((is>>std::ws).peek() != '}') {
    is >> key;
    is >> std::ws >> colon >> std::ws;
    is >> val;
    if ((is>>std::ws).peek() == ',') is >> comma;
    m[key] = val;
  }
  is >> std::ws >> paren;
  return is;
}


#endif /* end of include guard: PARAMETER_H */
