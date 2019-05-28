#ifndef RANDOM_H
#define RANDOM_H

class Rand
{ // Implementation of the highest quality recommended generator by "Numerical Recipes 3rd edition."
  // The constructor is called with an integer seed and creates an instance of the generator. The member functions
  // int64, doub, and int32 return the next values in the random sequence, as a variable type indicated by their names.
  // The period of the generator is ~ 3.138e57.
  //
  //  int64() = [Xorshift(LCG(u)) + Xorshift(v)] ^ MWC(w)
  //  LCG Modulo 2^64 for u.         (a = 2862933555777941757, c = 7046029254386353087)
  //  64-bit Xorshift for v.         (a1 = 17, a2 = 31, a3 = 8)
  //  Multiply with Carry with 2^32. (a = 4294957665)
  //  64-bit Xorshift for x from u.  (a1 = 21, a2 = 35, a3 = 4)
  private:
    unsigned long long  u, v, w;

  public:
    Rand(unsigned long long j);
    virtual ~Rand() { }

    unsigned long long  int64();
    double              doub();
    unsigned int        int32();
    int                 intIn(const int a, const int b);
};

Rand::Rand(unsigned long long j) : v(4101842887655102017LL), w(1) 
{ // Constructor. Call with any integer seed (except the value of v above). An example of a seed:
  //  unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  u = j^v;
  int64();
  v = u;
  int64();
  w = v;
  int64();
}

inline unsigned long long Rand::int64()
{ // Returns 64-bit random integer. 
  u = u * 2862933555777941757LL + 7046029254386353087LL;
  v ^= v >> 17;
  v ^= v << 31;
  v ^= v >> 8;
  w = 4294957665U*(w & 0xffffffff) + (w >> 32);
  unsigned long long x = u ^ (u << 21);
  x ^= x >> 35;
  x ^= x << 4;
  return (x + v) ^ w;
}

inline double Rand::doub()
{ // Returns random double-precision floating value in the range [0, 1]
  return 5.42101086242752217E-20 * int64();
}

inline unsigned int Rand::int32()
{ // Returns 32-bit random integer.
  return static_cast<unsigned int>(int64());
}

inline int Rand::intIn(const int a, const int b)
{ // Returns random int in the range [a,b).
  return static_cast<int>(a+int64()%(b-a));
}

#endif /* end of include guard: RANDOM_H */
