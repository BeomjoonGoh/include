#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <cmath>

namespace Util {
  class Logging
  { // Redirects clog to "logfileName".
    // Initiated by init() function, closed by, optional, call of terminate() function
    public:
      static std::ofstream logfileStream;
      static std::streambuf *old_rdbuf;
      static bool is_on;
    public:
      static void init(std::string logfileName);
      static void terminate();
  };
  std::ofstream   Logging::logfileStream;
  std::streambuf *Logging::old_rdbuf;
  bool            Logging::is_on = false;

  class Timer
  { // Timer class. Instantiate to start the timer and use elapsed() to get elapsed time.
    // Static function currentTD() records current time and date.
    private:
      std::chrono::time_point<std::chrono::system_clock> m_start;
    public:
      Timer() : m_start{std::chrono::system_clock::now()} { }
      virtual ~Timer() { }
      std::string elapsed();
      static std::string currentTD();
  };

  class outV
  { // Functor like class for ostream vector-like object v. Only instantiation is enough: cout << outV(v);
    private:
      int wid;
      int len;
      double *f;
    public:
      template <typename Container> outV(const Container &v, const int wid_=0, const char t='\0');
      virtual ~outV() { delete[] f; }
      friend std::ostream& operator<< (std::ostream &oS, const outV &A);
      std::string str();
  };

  std::ostream& operator<< (std::ostream &oS, const outV &A);
  void getComment(std::ifstream &inf, const std::string &inputf);
  void printHelp(const std::string &inputf);
  std::string ordinal(const int i);
  static const char ordinalSuffixes[][3] = {"th", "st", "nd", "rd"};

}

void Util::Logging::init(std::string logfileName)
{
  if (is_on) {
    std::cerr << "Util::Logging nested initialization." << std::endl;
    std::exit(1);
  }

  if (logfileName == "STDERR") return;
  is_on = true;
  logfileStream.open(logfileName);
  old_rdbuf = std::clog.rdbuf();          // Get the rdbuf of clog (We need it to reset the value before exiting).

  std::clog.rdbuf(logfileStream.rdbuf()); // Set the rdbuf of clog.
  std::clog << "Logging initiated.\n";
  std::atexit(Util::Logging::terminate);
}

void Util::Logging::terminate()
{
  if (!is_on) return;
  is_on = false;
  std::clog.flush();
  logfileStream.close();
  std::clog.rdbuf(old_rdbuf);             // Reset the rdbuf of clog.
}

std::string Util::Timer::elapsed()
{
  std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - m_start;

  double i0;
  double d = modf(elapsed_seconds.count(),&i0);
  int i = static_cast<int>(i0);

  int h = i/3600;
  int m = (i/60)%60;
  double s = i%60 + d;

  std::string out = "";
  if (h<10)
    out += "0";
  out += std::to_string(h)+":";
  if (m<10)
    out += "0";
  out += std::to_string(m)+":";
  if (s<10)
    out += "0";
  out += std::to_string(s);

  return out;
}

std::string Util::Timer::currentTD()
{
  auto now = std::chrono::system_clock::now();
  auto in_time_t = std::chrono::system_clock::to_time_t(now);

  char buf[100];
  std::strftime(buf, sizeof(buf), "%Y-%m-%d %X", std::localtime(&in_time_t));
  //std::stringstream ss;
  //ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
  //return ss.str();
  return buf;
}

template <typename Container>
inline Util::outV::outV(const Container &v, const int wid_, const char t)
{
  wid = wid_;
  len = v.size();
  f = new double[len];
  for (int i = 0; i < len; i++)
    f[i] = (t=='i') ? static_cast<int>(v[i]) : v[i];
}

inline std::ostream& Util::operator<< (std::ostream &oS, const outV &A)
{
  for (int i = 0; i < A.len-1; i++)
    oS << std::setw(A.wid) << A.f[i] << " ";
  oS << std::setw(A.wid) << A.f[A.len-1];
  return oS;
}

std::string Util::outV::str()
{
  std::ostringstream oss;
  oss << *this;
  return oss.str();
}

void Util::getComment(std::ifstream &inf, const std::string &inputf)
{
  inf.clear();                    // forget we hit the end of file
  inf.seekg(0, std::ios::beg);    // move to the start of the file
  if (inf.peek() == '#') {
    std::clog << "Comment from " << inputf << " :\n";
    std::string strInput;
    getline(inf, strInput);
    std::clog << "   " << strInput << '\n';
  }
}

void Util::printHelp(const std::string &inputf)
{
  std::ifstream inf{inputf};
  if (inf.good()) std::cerr << inf.rdbuf();
  else            std::cerr << "Cannot find the help document." << std::endl;
}

std::string Util::ordinal(const int i)
{
  int ord = i % 100;
  if (ord / 10 == 1)
    ord = 0;
  ord %= 10;
  if (ord > 3)
    ord = 0;
  return std::to_string(i)+std::string{ordinalSuffixes[ord]};
}


#endif /* end of include guard: UTIL_H */
