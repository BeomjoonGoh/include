#ifndef PARAM_H
#define PARAM_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "util.h"

using namespace std;

class Param
{ // To add a new parameter,
  //    1) Define public static member variable,
  //    2) Externalize it,
  //    3) Set its default value in setDefaults(),
  //    4) Make it possible to be read from inputf,     (optional)
  //    5) Let it be printed to log.                    (optional)      
public:
  static double       beta;
  static double       U;
  static int          Niom;
  static string  logfile;

  static void ReadParams(string inputf);
  static void Log();

private:
  static bool m_exist;
  static string m_inputf;

  static void setDefaults()
  {
    beta        =     300.0;
    U           =       3.0;
    Niom        =    2000;
    logfile     = "log.txt";
  }
};
extern double       Param::beta;
extern double       Param::U;
extern int          Param::Niom;
extern string  Param::logfile;

extern bool Param::m_exist;
extern string Param::m_inputf;

void Param::ReadParams(string inputf)
{
  setDefaults();

  ifstream inf {inputf};
  if (inf.good()) {
    m_exist = true;
    m_inputf = inputf;
    while (inf) {
      string strInput;
      inf >> strInput;
      if (strInput[0]=='#') inf.ignore(2000,'\n');
      if (strInput=="beta")                           inf >> beta;
      if (strInput=="U")                              inf >> U;
      if (strInput=="Niom")                           inf >> Niom;
      if (strInput=="logfile")                        inf >> logfile;
    }
  } else {
    m_exist = false;
    cerr << RED("Missing input file!") << "Proceeds only with default values." << endl;
  }
}

void Param::Log()
{
  clog << "============================================================\n";
  clog << "Input parameters:\n";
  if (m_exist) clog << "    Read from " << m_inputf << " file, set defaults otherwise.\n";
  else         clog << "    Missing input file! Proceeds only with default values.\n";

  clog << "  beta          = " << beta         << '\n'
       << "  U             = " << U            << '\n'
       << "  Niom          = " << Niom         << '\n'
       << "============================================================\n";
}


#endif /* end of include guard: PARAM_H */
