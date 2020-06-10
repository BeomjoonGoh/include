#ifndef ASSERT_H
#define ASSERT_H

#include <iostream>
#include <cstdlib>

#define RED(mess)       "\033[01;38;5;210m" << (mess) << " \033[0m"
#define GREEN(mess)     "\033[01;38;5;114m" << (mess) << " \033[0m"
#define YELLOW(mess)    "\033[01;38;5;186m" << (mess) << " \033[0m"
#define BLUE(mess)      "\033[01;38;5;117m" << (mess) << " \033[0m"
#define PURPLE(mess)    "\033[01;38;5;211m" << (mess) << " \033[0m"
#define WHITE(mess)     "\033[01;38;5;231m" << (mess) << " \033[0m"
#define DEFAULTCOLOR    "\033[0m"

#define quit() {\
  std::cerr << BLUE("QUIT HERE >")<<__FILE__<<":"<<__LINE__<<" In function '"<<__func__<<"'"<< std::endl;\
  std::exit(EXIT_FAILURE);\
}

#define quitif(cond, mess) {\
  if ((cond)) { \
    std::cerr << RED("QUIT!") << "since " << "(" #cond ")\n" \
    << __FILE__ << ":" << __LINE__ << " In function '" << __func__ << "': " << mess << std::endl; \
    std::exit(EXIT_FAILURE);\
  }\
}

#ifndef NDEBUG
  #define assert(cond, mess) {\
    if (!(cond)) { \
      std::cerr << RED("Assertion!") << "(" #cond ") failed. \n" \
                << __FILE__ << ":" << __LINE__ << " In function '" << __func__ << "': " << mess << std::endl; \
      std::exit(EXIT_FAILURE);\
    }\
  }

  #define assert0(cond) {\
    if (!(cond)) { \
      std::cerr << RED("Assertion!") << "(" #cond ") failed. \n" \
                << __FILE__ << ":" << __LINE__ << " In function '" << __func__ << "'" << std::endl; \
      std::exit(EXIT_FAILURE);\
    }\
  }
#else
  #define assert(cond, mess)
  #define assert0(cond)
#endif // NDEBUG

#endif /* end of include guard: ASSERT_H */
