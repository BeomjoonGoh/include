#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <cmath>
#include <unistd.h>
#include <sys/wait.h>

#include "assert.h"

namespace Util {
  class Logging
  { // Redirects clog to "logfileName".
    // Initiated by init() function, closed by, optional, call of terminate() function
    public:
      static std::ofstream logfileStream;
      static std::streambuf *old_rdbuf;
      static bool is_on;
    public:
      static void init(std::string logfileName, std::ios_base::openmode mode = std::ios_base::out);
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
      ~Timer() { }
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
      ~outV() { delete[] f; }
      friend std::ostream& operator<< (std::ostream &oS, const outV &A);
      std::string str();
  };
  std::ostream& operator<< (std::ostream &oS, const outV &A);

  void runCommand(const std::string &command);

  void getComment(std::ifstream &inf, const std::string &inputf);
  std::string unindent(const char *p);
  std::string ordinal(const int i);
  int wordCount(const std::string &text);
}

void Util::Logging::init(std::string logfileName, std::ios_base::openmode mode)
{
  quitif(is_on, "Util::Logging nested initialization.");

  if (logfileName == "STDERR") return;
  is_on = true;
  logfileStream.open(logfileName, mode);
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
  double d = std::modf(elapsed_seconds.count(),&i0);
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

void Util::runCommand(const std::string &command)
{
  std::clog << "Util::runCommand(\"" << command << "\")\n";
  std::cout.flush();
  std::clog.flush();

  int size = Util::wordCount(command);
  std::string pieces[size];
  const char* argv[size+1];

  std::stringstream s(command);
  for (int i = 0; s >> std::ws >> pieces[i]; i++)
    argv[i] = pieces[i].c_str();
  argv[size] = NULL;

  pid_t pid = fork();
  if (pid == 0) { // child
    execv(argv[0], const_cast<char* const *>(argv));
    quitif(true, RED("child") << ": execv failed to run");
  } else if (pid > 0) { // parent 
    int status_child;
    if (waitpid(pid, &status_child, 0) > 0) {
      if (WIFEXITED(status_child)) {
        switch (WEXITSTATUS(status_child)) { 
          case 0:   { std::clog << "parent: process terminated with zero returned." << std::endl; break; }
          case 127: quitif(true, RED("parent") << ": execv failed.");
          default:  quitif(true, RED("parent") << ": process terminated with non-zero returned.");
        }
      } else
        quitif(true, RED("parent") << ": process terminated abnormally.");
    } else
      quitif(true, RED("parent") << ": waitpid failed.");
  } else
    quitif(true, "Failed to fork.");
}

void Util::getComment(std::ifstream &inf, const std::string &inputf)
{
  inf.clear();                    // forget we hit the end of file
  inf.seekg(0, std::ios::beg);    // move to the start of the file
  if (inf.peek() == '#') {
    std::clog << "Comment from " << inputf << " :\n";
    std::string strInput;
    std::getline(inf, strInput);
    std::clog << "   " << strInput << '\n';
  }
}

std::string Util::unindent(const char *p)
{
  if (*p == '\n') p++;

  std::string indents;
  for (; std::isspace(*p) && *p != '\n'; p++)
    indents += *p;

  std::string result;
  for (; *p; p++) {
    result += *p;
    if (*p == '\n' && std::string(p+1, 0, indents.size()) == indents)
      p += indents.size();
  }
  return result;
}

std::string Util::ordinal(const int i)
{
  const char ordinalSuffixes[][3] = {"th", "st", "nd", "rd"};
  int ord = i % 100;
  if (ord / 10 == 1)
    ord = 0;
  ord %= 10;
  if (ord > 3)
    ord = 0;
  return std::to_string(i)+std::string{ordinalSuffixes[ord]};
}

int Util::wordCount(const std::string &text)
{
  int size = 0;
  std::string word;
  for (std::stringstream s(text); s >> std::ws >> word; size++);
  return size;
}


#endif /* end of include guard: UTIL_H */
