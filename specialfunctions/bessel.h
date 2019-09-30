#ifndef BESSEL_H
#define BESSEL_H

#include "maths.h"

class BesselN
{ // Bessel function of integer order n of the first  kind Jn(n,x), n >= 0, any x.
  // Bessel function of integer order n of the second kind Yn(n,x), n >= 0, x > 0.
  // Modified Bessel function of integer order n of the first  kind In(n,x), n >= 0, any x.
  // Modified Bessel function of integer order n of the second kind Kn(n,x), n >= 0, x > 0.
  public:
    double Jn(const int n, const double x);
    double Yn(const int n, const double x);

    double In(const int n, const double x);
    double Kn(const int n, const double x);

  private:
    double J0(const double x);
    double J1(const double x);
    double Y0(const double x);
    double Y1(const double x);
    void rational(const double x, const double *r, const double *s, const int n);
    void asymptot(const double *pn, const double *pd, const double *qn, const double *qd, const double fac);

    double I0(const double x);
    double I1(const double x);
    double K0(const double x);
    double K1(const double x);
    double polynomial(const double *cof, const int n, const double x);

  private:
    static const double xJ00, xJ10, xJ01, xJ11, twoopi, pio4;
    static const double J0r[7], J0s[7], J0pn[5], J0pd[5], J0qn[5], J0qd[5];
    static const double J1r[7], J1s[7], J1pn[5], J1pd[5], J1qn[5], J1qd[5];
    static const double Y0r[9], Y0s[9], Y0pn[5], Y0pd[5], Y0qn[5], Y0qd[5];
    static const double Y1r[8], Y1s[8], Y1pn[5], Y1pd[5], Y1qn[5], Y1qd[5];

    static const double I0p[14], I0q[5], I0pp[5], I0qq[6];
    static const double I1p[14], I1q[5], I1pp[5], I1qq[6];
    static const double K0pi[5], K0qi[3], K0p[5], K0q[3], K0pp[8], K0qq[8];
    static const double K1pi[5], K1qi[3], K1p[5], K1q[3], K1pp[8], K1qq[8];

    const int IEXP = std::numeric_limits<double>::max_exponent/2;
    double nump, denp, numq, denq, y, z, ax, xx;
};

class Bessel : public BesselN
{ // Object for Bessel functions of arbitrary order nu and related functions.
  public:
    Bessel() : xJY(9.99e99), nuJY(9.99e99), xIK(9.99e99), nuIK(9.99e99), xAB(9.99e99), xjy(9.99e99), njy(-9999) { }
    virtual ~Bessel() { }

    double Jnu(const double nu, const double x);
    double Ynu(const double nu, const double x);
    double Inu(const double nu, const double x);
    double Knu(const double nu, const double x);
    double Ai(const double x);
    double Bi(const double x);
    double jn(const int n, const double x);
    double yn(const int n, const double x);

  private:
    void besselJY(const double nu, const double x);
    void besselIK(const double nu, const double x);
    void airy(const double x);
    void spherical(const int n, const double x);

    inline double chebevshev(const double *c, const int m, const double x)
    { // Utility used by besselJY and besselIK, evaluates Chebyshev series.
      double d = 0.0, dd = 0.0;
      for (int j = m-1; j > 0; j--) {
        double sv = d;
        d = 2.*x*d-dd+c[j];
        dd = sv;
      }
      return x*d-dd+0.5*c[0];
    }

  private:
    double Jo, Yo, Jpo, Ypo;
    double Io, Ko, Ipo, Kpo;
    double Aio,Bio,Aipo,Bipo;
    double jo, yo, jpo, ypo;

    double xJY, nuJY, xIK, nuIK, xAB, xjy;
    int njy;


    static const int MAXIT = 10000;
    static const int NUSE1 = 7, NUSE2 = 8;
    static const double c1[NUSE1], c2[NUSE2];
    const double Fmin = std::numeric_limits<double>::min()/Maths::absEpsilon;
    const double XMIN = 2.0;

};

double BesselN::J0(const double x)
{
  ax = Maths::abs(x);
  if (ax < 8.0) {
    rational(x, J0r, J0s, 6);
    return nump*(y-xJ00)*(y-xJ10)/denp;
  } else {
    asymptot(J0pn, J0pd, J0qn, J0qd, 1.);
    return std::sqrt(twoopi/ax)*(std::cos(xx)*nump/denp-z*std::sin(xx)*numq/denq);
  }
}

double BesselN::J1(const double x)
{
  ax = Maths::abs(x);
  if (ax < 8.0) {
    rational(x, J1r, J1s, 6);
    return x*nump*(y-xJ01)*(y-xJ11)/denp;
  } else {
    asymptot(J1pn, J1pd, J1qn, J1qd, 3.);
    double ans = std::sqrt(twoopi/ax)*(std::cos(xx)*nump/denp-z*std::sin(xx)*numq/denq);
    return x > 0.0 ? ans : -ans;
  }
}

double BesselN::Y0(const double x)
{
  if (x < 8.0) {
    double J0x = J0(x);
    rational(x, Y0r, Y0s, 8);
    return nump/denp+twoopi*J0x*std::log(x);
  } else {
    ax = x;
    asymptot(Y0pn, Y0pd, Y0qn, Y0qd, 1.);
    return std::sqrt(twoopi/x)*(std::sin(xx)*nump/denp+z*std::cos(xx)*numq/denq);
  }
}

double BesselN::Y1(const double x)
{
  if (x < 8.0) {
    double J1x = J1(x);
    rational(x, Y1r, Y1s, 7);
    return x*nump/denp+twoopi*(J1x*std::log(x)-1.0/x);
  } else {
    ax = x;
    asymptot(Y1pn, Y1pd, Y1qn, Y1qd, 3.);
    return std::sqrt(twoopi/x)*(std::sin(xx)*nump/denp+z*std::cos(xx)*numq/denq);
  }
}

double BesselN::Yn(const int n, const double x)
{
  if (n == 0) return Y0(x);
  if (n == 1) return Y1(x);
  double by, bym, byp;
  by = Y1(x);
  bym = Y0(x);
  double tox = 2.0/x;
  for (int j = 1; j < n; j++) { // Recurrence.
    byp = j*tox*by-bym;
    bym = by;
    by = byp;
  }
  return by;
}

double BesselN::Jn(const int n, const double x)
{
  if (n == 0) return J0(x);
  if (n == 1) return J1(x);

  const double ACC = 160.0; // ACC determines accuracy.
  double ax = Maths::abs(x);

  double bj, bjm, bjp;
  double ans;
  if (ax*ax <= 8.0*Maths::doubleMin) {
    return 0.0;
  } else if (ax > double(n)) { // Upwards recurrence from J0 and J1.
    double tox = 2.0/ax;
    bjm = J0(ax);
    bj = J1(ax);
    for (int j = 1; j < n; j++) {
      bjp = j*tox*bj-bjm;
      bjm = bj;
      bj = bjp;
    }
    ans = bj;
  } else { // Downward recurrence from an even m here computed.
    double tox = 2.0/ax;
    // jsum will alternate between false and true; when it is true, we accumulate in sum the even terms in (5.4.16).
    bool jsum = false;
    bjp = ans = 0.0;
    bj = 1.0;
    int k;
    double sum = 0.0;
    int m = (n+static_cast<int>(std::sqrt(ACC*n))) >> 1 << 1;
    for (int j = m; j > 0; j--) { // The downward recurrence.
      bjm = j*tox*bj-bjp;
      bjp = bj;
      bj = bjm;
      std::frexp(bj, &k);
      if (k > IEXP) { // Renormalize to prevent overflows.
        bj =  std::ldexp(bj, -IEXP);
        bjp = std::ldexp(bjp, -IEXP);
        ans = std::ldexp(ans, -IEXP);
        sum = std::ldexp(sum, -IEXP);
      }
      if (jsum) sum += bj; // Accumulate the sum.
      jsum = !jsum; // Change false to true or vice versa.
      if (j == n) ans = bjp; // Save the unnormalized answer.
    }
    // Compute (5.4.16) and use it to normalize the answer.
    sum = 2.0*sum-bj;
    ans /= sum;
  }

  return x < 0.0 && (n & 1) ? -ans : ans;
}

inline void BesselN::rational(const double x, const double *r, const double *s, const int n)
{ // Evaluates rational approximation.
  y = x*x;
  z = 64.0-y;
  nump = r[n];
  denp = s[n];
  for (int i = n-1; i >= 0; i--) {
    nump = nump*z+r[i];
    denp = denp*y+s[i];
  }
}

inline void BesselN::asymptot(const double *pn, const double *pd, const double *qn, const double *qd, const double fac)
{ // Evaluates asymptotic approximation.
  z = 8.0/ax;
  y = z*z;
  xx = ax-fac*pio4;
  nump = pn[4];
  denp = pd[4];
  numq = qn[4];
  denq = qd[4];
  for (int i = 3; i >= 0; i--) {
    nump = nump*y+pn[i];
    denp = denp*y+pd[i];
    numq = numq*y+qn[i];
    denq = denq*y+qd[i];
  }
}

double BesselN::I0(const double x)
{
  ax = Maths::abs(x);
  if (ax < 15.0) {
    y = x*x;
    return polynomial(I0p, 13, y)/polynomial(I0q, 4, 225.-y);
  } else {
    z = 1.0-15.0/ax;
    return std::exp(ax)*polynomial(I0pp, 4, z)/(polynomial(I0qq, 5, z)*std::sqrt(ax));
  }
}

double BesselN::I1(const double x)
{
  ax = Maths::abs(x);
  if (ax < 15.0) {
    y = x*x;
    return x*polynomial(I1p, 13, y)/polynomial(I1q, 4, 225.-y);
  } else {
    z = 1.0-15.0/ax;
    double ans = std::exp(ax)*polynomial(I1pp, 4, z)/(polynomial(I1qq, 5, z)*std::sqrt(ax));
    return x > 0.0 ? ans : -ans;
  }
}

double BesselN::K0(const double x)
{
  if (x <= 1.0) {
    z = x*x;
    double term = polynomial(K0pi, 4, z)*std::log(x)/polynomial(K0qi, 2, 1.-z);
    return polynomial(K0p, 4, z)/polynomial(K0q, 2, 1.-z)-term;
  } else {
    z = 1.0/x;
    return std::exp(-x)*polynomial(K0pp, 7, z)/(polynomial(K0qq, 7, z)*std::sqrt(x));
  }
}

double BesselN::K1(const double x)
{
  if (x <= 1.0) {
    z = x*x;
    double term = polynomial(K1pi, 4, z)*std::log(x)/polynomial(K1qi, 2, 1.-z);
    return x*(polynomial(K1p, 4, z)/polynomial(K1q, 2, 1.-z)+term)+1./x;
  } else {
    z = 1.0/x;
    return std::exp(-x)*polynomial(K1pp, 7, z)/(polynomial(K1qq, 7, z)*std::sqrt(x));
  }
}

double BesselN::Kn(const int n, const double x)
{
  if (n == 0) return K0(x);
  if (n == 1) return K1(x);
  double bk, bkm, bkp;
  bkm = K0(x);
  bk = K1(x);
  double tox = 2.0/x;
  for (int j = 1; j < n; j++) {
    bkp = bkm+j*tox*bk;
    bkm = bk;
    bk = bkp;
  }
  return bk;
}

double BesselN::In(const int n, const double x)
{
  if (n == 0) return I0(x);
  if (n == 1) return I1(x);

  const double ACC = 200.0; // ACC determines accuracy.

  double bi, bim, bip;
  double ans;
  if (x*x <= 8.0*Maths::doubleMin) {
    return 0.0;
  } else {
    double tox = 2.0/Maths::abs(x);
    bip = ans = 0.0;
    bi = 1.0;
    int k;
    int m = (n+static_cast<int>(std::sqrt(ACC*n))) << 1;
    for (int j = m; j > 0; j--) { // Downward recurrence.
      bim = bip+j*tox*bi;
      bip = bi;
      bi = bim;
      std::frexp(bi, &k);
      if (k > IEXP) { // Renormalize to prevent overflows.
        ans = std::ldexp(ans, -IEXP);
        bi  = std::ldexp(bi, -IEXP);
        bip = std::ldexp(bip, -IEXP);
      }
      if (j == n) ans = bip;
    }
    ans *= I0(x)/bi; // Normalize with bessi0.
    return x < 0.0 && (n & 1) ? -ans : ans;
  }
}

inline double BesselN::polynomial(const double *cof, const int n, const double x)
{ // Evaluate a polynomial.
  double ans = cof[n];
  for (int i = n-1; i >= 0; i--)
    ans = ans*x+cof[i];
  return ans;
}

// First two zeros of J0, 2/pi and pi/4
const double BesselN::xJ00 = 5.783185962946785;
const double BesselN::xJ10 = 3.047126234366209e1;
// First two zeros of J1, 2/pi and pi/4
const double BesselN::xJ01 = 1.468197064212389e1;
const double BesselN::xJ11 = 4.921845632169460e1;
const double BesselN::twoopi = 0.6366197723675813;
const double BesselN::pio4 = 0.7853981633974483;

const double BesselN::J0r[]  = {1.682397144220462e-4, 2.058861258868952e-5, 5.288947320067750e-7, 5.557173907680151e-9, 2.865540042042604e-11, 7.398972674152181e-14, 7.925088479679688e-17};
const double BesselN::J0s[]  = {1.0, 1.019685405805929e-2, 5.130296867064666e-5, 1.659702063950243e-7, 3.728997574317067e-10, 5.709292619977798e-13, 4.932979170744996e-16};
const double BesselN::J0pn[] = {9.999999999999999e-1, 1.039698629715637, 2.576910172633398e-1, 1.504152485749669e-2, 1.052598413585270e-4};
const double BesselN::J0pd[] = {1.0, 1.040797262528109, 2.588070904043728e-1, 1.529954477721284e-2, 1.168931211650012e-4};
const double BesselN::J0qn[] = {-1.562499999999992e-2, -1.920039317065641e-2, -5.827951791963418e-3, -4.372674978482726e-4, -3.895839560412374e-6};
const double BesselN::J0qd[] = {1.0, 1.237980436358390, 3.838793938147116e-1, 3.100323481550864e-2, 4.165515825072393e-4};
const double BesselN::J1r[]  = {7.309637831891357e-5, 3.551248884746503e-6, 5.820673901730427e-8, 4.500650342170622e-10, 1.831596352149641e-12, 3.891583573305035e-15, 3.524978592527982e-18};
const double BesselN::J1s[]  = {1.0, 9.398354768446072e-3, 4.328946737100230e-5, 1.271526296341915e-7, 2.566305357932989e-10, 3.477378203574266e-13, 2.593535427519985e-16};
const double BesselN::J1pn[] = {1.0, 1.014039111045313, 2.426762348629863e-1, 1.350308200342000e-2, 9.516522033988099e-5};
const double BesselN::J1pd[] = {1.0, 1.012208056357845, 2.408580305488938e-1, 1.309511056184273e-2, 7.746422941504713e-5};
const double BesselN::J1qn[] = {4.687499999999991e-2, 5.652407388406023e-2, 1.676531273460512e-2, 1.231216817715814e-3, 1.178364381441801e-5};
const double BesselN::J1qd[] = {1.0, 1.210119370463693, 3.626494789275638e-1, 2.761695824829316e-2, 3.240517192670181e-4};
const double BesselN::Y0r[]  = {-7.653778457189104e-3, -5.854760129990403e-2, 3.720671300654721e-4, 3.313722284628089e-5, 4.247761237036536e-8, -4.134562661019613e-9, -3.382190331837473e-11, -1.017764126587862e-13, -1.107646382675456e-16};
const double BesselN::Y0s[]  = {1.0, 1.125494540257841e-2, 6.427210537081400e-5, 2.462520624294959e-7, 7.029372432344291e-10, 1.560784108184928e-12, 2.702374957564761e-15, 3.468496737915257e-18, 2.716600180811817e-21};
const double BesselN::Y0pn[] = {9.999999999999999e-1, 1.039698629715637, 2.576910172633398e-1, 1.504152485749669e-2, 1.052598413585270e-4};
const double BesselN::Y0pd[] = {1.0, 1.040797262528109, 2.588070904043728e-1, 1.529954477721284e-2, 1.168931211650012e-4};
const double BesselN::Y0qn[] = {-1.562499999999992e-2, -1.920039317065641e-2, -5.827951791963418e-3, -4.372674978482726e-4, -3.895839560412374e-6};
const double BesselN::Y0qd[] = {1.0, 1.237980436358390, 3.838793938147116e-1, 3.100323481550864e-2, 4.165515825072393e-4};
const double BesselN::Y1r[]  = {-1.041835425863234e-1, -1.135093963908952e-5, 2.212118520638132e-4, 1.270981874287763e-6, -3.982892100836748e-8, -4.820712110115943e-10, -1.929392690596969e-12, -2.725259514545605e-15};
const double BesselN::Y1s[]  = {1.0, 1.186694184425838e-2, 7.121205411175519e-5, 2.847142454085055e-7, 8.364240962784899e-10, 1.858128283833724e-12, 3.018846060781846e-15, 3.015798735815980e-18};
const double BesselN::Y1pn[] = {1.0, 1.014039111045313, 2.426762348629863e-1, 1.350308200342000e-2, 9.516522033988099e-5};
const double BesselN::Y1pd[] = {1.0, 1.012208056357845, 2.408580305488938e-1, 1.309511056184273e-2, 7.746422941504713e-5};
const double BesselN::Y1qn[] = {4.687499999999991e-2, 5.652407388406023e-2, 1.676531273460512e-2, 1.231216817715814e-3, 1.178364381441801e-5};
const double BesselN::Y1qd[] = {1.0, 1.210119370463693, 3.626494789275638e-1, 2.761695824829316e-2, 3.240517192670181e-4};

const double BesselN::I0p[]  = {9.999999999999997e-1, 2.466405579426905e-1, 1.478980363444585e-2, 3.826993559940360e-4, 5.395676869878828e-6, 4.700912200921704e-8, 2.733894920915608e-10, 1.115830108455192e-12, 3.301093025084127e-15, 7.209167098020555e-18, 1.166898488777214e-20, 1.378948246502109e-23, 1.124884061857506e-26, 5.498556929587117e-30};
const double BesselN::I0q[]  = {4.463598170691436e-1, 1.702205745042606e-3, 2.792125684538934e-6, 2.369902034785866e-9, 8.965900179621208e-13};
const double BesselN::I0pp[] = {1.192273748120670e-1, 1.947452015979746e-1, 7.629241821600588e-2, 8.474903580801549e-3, 2.023821945835647e-4};
const double BesselN::I0qq[] = {2.962898424533095e-1, 4.866115913196384e-1, 1.938352806477617e-1, 2.261671093400046e-2, 6.450448095075585e-4, 1.529835782400450e-6};
const double BesselN::I1p[]  = {5.000000000000000e-1, 6.090824836578078e-2, 2.407288574545340e-3, 4.622311145544158e-5, 5.161743818147913e-7, 3.712362374847555e-9, 1.833983433811517e-11, 6.493125133990706e-14, 1.693074927497696e-16, 3.299609473102338e-19, 4.813071975603122e-22, 5.164275442089090e-25, 3.846870021788629e-28, 1.712948291408736e-31};
const double BesselN::I1q[]  = {4.665973211630446e-1, 1.677754477613006e-3, 2.583049634689725e-6, 2.045930934253556e-9, 7.166133240195285e-13};
const double BesselN::I1pp[] = {1.286515211317124e-1, 1.930915272916783e-1, 6.965689298161343e-2, 7.345978783504595e-3, 1.963602129240502e-4};
const double BesselN::I1qq[] = {3.309385098860755e-1, 4.878218424097628e-1, 1.663088501568696e-1, 1.473541892809522e-2, 1.964131438571051e-4, -1.034524660214173e-6};
const double BesselN::K0pi[] = {1.0, 2.346487949187396e-1, 1.187082088663404e-2, 2.150707366040937e-4, 1.425433617130587e-6};
const double BesselN::K0qi[] = {9.847324170755358e-1, 1.518396076767770e-2, 8.362215678646257e-5};
const double BesselN::K0p[]  = {1.159315156584126e-1, 2.770731240515333e-1, 2.066458134619875e-2, 4.574734709978264e-4, 3.454715527986737e-6};
const double BesselN::K0q[]  = {9.836249671709183e-1, 1.627693622304549e-2, 9.809660603621949e-5};
const double BesselN::K0pp[] = {1.253314137315499, 1.475731032429900e1, 6.123767403223466e1, 1.121012633939949e2, 9.285288485892228e1, 3.198289277679660e1, 3.595376024148513, 6.160228690102976e-2};
const double BesselN::K0qq[] = {1.0, 1.189963006673403e1, 5.027773590829784e1, 9.496513373427093e1, 8.318077493230258e1, 3.181399777449301e1, 4.443672926432041, 1.408295601966600e-1};
const double BesselN::K1pi[] = {0.5, 5.598072040178741e-2, 1.818666382168295e-3, 2.397509908859959e-5, 1.239567816344855e-7};
const double BesselN::K1qi[] = {9.870202601341150e-1, 1.292092053534579e-2, 5.881933053917096e-5};
const double BesselN::K1p[]  = {-3.079657578292062e-1, -8.109417631822442e-2, -3.477550948593604e-3, -5.385594871975406e-5, -3.110372465429008e-7};
const double BesselN::K1q[]  = {9.861813171751389e-1, 1.375094061153160e-2, 6.774221332947002e-5};
const double BesselN::K1pp[] = {1.253314137315502, 1.457171340220454e1, 6.063161173098803e1, 1.147386690867892e2, 1.040442011439181e2, 4.356596656837691e1, 7.265230396353690, 3.144418558991021e-1};
const double BesselN::K1qq[] = {1.0, 1.125154514806458e1, 4.427488496597630e1, 7.616113213117645e1, 5.863377227890893e1, 1.850303673841586e1, 1.857244676566022, 2.538540887654872e-2};



double Bessel::Jnu(const double nu, const double x)
{
  if (nu != nuJY || x != xJY) besselJY(nu, x);
  return Jo;
}

double Bessel::Ynu(const double nu, const double x)
{
  if (nu != nuJY || x != xJY) besselJY(nu, x);
  return Yo;
}

double Bessel::Inu(const double nu, const double x)
{
  if (nu != nuIK || x != xIK) besselIK(nu, x);
  return Io;
}

double Bessel::Knu(const double nu, const double x)
{
  if (nu != nuIK || x != xIK) besselIK(nu, x);
  return Ko;
}

double Bessel::Ai(const double x)
{
  if (x != xAB) airy(x);
  return Aio;
}

double Bessel::Bi(const double x)
{
  if (x != xAB) airy(x);
  return Bio;
}

double Bessel::jn(const int n, const double x)
{
  if (n != njy || x != xjy) spherical(n, x);
  return jo;
}

double Bessel::yn(const int n, const double x)
{
  if (n != njy || x != xjy) spherical(n, x);
  return yo;
}

void Bessel::besselJY(const double nu, const double x)
{ // Sets Jo, Yo, Jpo, and Ypo respectively to the Bessel functions Jnu(x), Ynu(x) and their derivatives Jnu'(x),
  // Ynu'(x), for positive x and for xnu = \nu >= 0. The relative accuracy is within one or two significant digits of
  // Maths::absEpsilon, except near a zero of one of the functions, where Maths::absEpsilon controls its absolute accuracy. Fmin is a number close
  // to the machine's smallest floating-point number.
  double rjmu, ry1, rymu;
  int i;
  if (x <= 0.0 || nu < 0.0) throw("bad arguments in besselJY");
  // nl is the number of downward recurrences of the J's and upward recurrences of Y's. xmu lies between -1/2 and 1/2
  // for x < XMIN, while it is chosen so that x is greater than the turning point for x >= XMIN.
  int nl = (x < XMIN ? static_cast<int>(nu+0.5) : Maths::max(static_cast<int>(nu-x+1.5), 0));
  double xmu = nu-nl;
  double xmu2 = xmu*xmu;
  double xi = 1.0/x;
  double xi2 = 2.0*xi;
  double w = xi2/Maths::pi; // The Wronskian. Evaluate CF1 by modified Lentz's method (÷5.2)
  int isign = 1; // isign keeps track of sign changes in the denominator.
  double h = nu*xi;
  if (h < Fmin) h = Fmin;
  double b = xi2*nu;
  double d = 0.0;
  double c = h;
  for (i = 0; i < MAXIT; i++) {
    b += xi2;
    d = b-d;
    if (Maths::abs(d) < Fmin) d = Fmin;
    c = b-1.0/c;
    if (Maths::abs(c) < Fmin) c = Fmin;
    d = 1.0/d;
    double del = c*d;
    h = del*h;
    if (d < 0.0) isign = -isign;
    if (Maths::abs(del-1.0) <= Maths::absEpsilon) break;
  }
  if (i >= MAXIT)
    throw("x too large in besselJY; try asymptotic expansion");
  double rjl = isign*Fmin; // Initialize Jnu and Jnu0 for downward recurrence.
  double rjpl = h*rjl;
  double rjl1 = rjl; // Store values for later rescaling.
  double rjp1 = rjpl;
  double fact = nu*xi;
  for (int l = nl-1; l >= 0; l--) {
    double rjtemp = fact*rjl+rjpl;
    fact -= xi;
    rjpl = fact*rjtemp-rjl;
    rjl = rjtemp;
  }
  if (rjl == 0.0) rjl = Maths::absEpsilon;
  double f = rjpl/rjl; // Now have unnormalized Jmu and Jmu' .

  if (x < XMIN) { // Use series.
    double x2 = 0.5*x;
    double pimu = Maths::pi*xmu;
    fact = (Maths::abs(pimu) < Maths::absEpsilon ? 1.0 : pimu/std::sin(pimu));
    d = -std::log(x2);
    double e = xmu*d;
    double fact2 = (Maths::abs(e) < Maths::absEpsilon ? 1.0 : std::sinh(e)/e);
    double xx = 8.0*xmu*xmu-1.0; // Evaluates Gamma_1 and Gamma_2 by Chebyshev expansion for |x| <= 1/2. Also returns 1/Gamma(1 \pm x).
    double gam1 = chebevshev(c1, NUSE1, xx);
    double gam2 = chebevshev(c2, NUSE2, xx);
    double ff = 2.0/Maths::pi*fact*(gam1*std::cosh(e)+gam2*fact2*d); // f0.
    e = std::exp(e);
    double p = e/((gam2-xmu*gam1)*Maths::pi); // p0.
    double q = 1.0/(e*Maths::pi*(gam2+xmu*gam1)); // q0.
    double pimu2 = 0.5*pimu;
    double fact3 = (Maths::abs(pimu2) < Maths::absEpsilon ? 1.0 : std::sin(pimu2)/pimu2);
    double r = Maths::pi*pimu2*fact3*fact3;
    c = 1.0;
    d = -x2*x2;
    double sum = ff+r*q;
    double sum1 = p;
    for (i = 1; i <= MAXIT; i++) {
      ff = (i*ff+p+q)/(i*i-xmu2);
      c *= (d/i);
      p /= (i-xmu);
      q /= (i+xmu);
      double del = c*(ff+r*q);
      sum += del;
      sum1 += c*p-i*del;
      if (Maths::abs(del) < (1.0+Maths::abs(sum))*Maths::absEpsilon) break;
    }
    if (i > MAXIT) throw("bessy series failed to converge");
    rymu = -sum;
    ry1 = -sum1*xi2;
    double rymup = xmu*xi*rymu-ry1; // Equation (6.6.13)
    rjmu = w/(rymup-f*rymu);
  } else { // EvaluateCF2bymodifiedLentz'smethod(÷5.2).
    double a = 0.25-xmu2;
    double p = -0.5*xi;
    double q = 1.0;
    double br = 2.0*x;
    double bi = 2.0;
    fact = a*xi/(p*p+q*q);
    double cr = br+q*fact;
    double ci = bi+p*fact;
    double den = br*br+bi*bi;
    double dr = br/den;
    double di = -bi/den;
    double dlr = cr*dr-ci*di;
    double dli = cr*di+ci*dr;
    double temp = p*dlr-q*dli;
    q = p*dli+q*dlr;
    p = temp;
    for (i = 1; i < MAXIT; i++) {
      a += 2*i;
      bi += 2.0;
      dr = a*dr+br;
      di = a*di+bi;
      if (Maths::abs(dr)+Maths::abs(di) < Fmin) dr = Fmin;
      fact = a/(cr*cr+ci*ci);
      cr = br+cr*fact;
      ci = bi-ci*fact;
      if (Maths::abs(cr)+Maths::abs(ci) < Fmin) cr = Fmin;
      den = dr*dr+di*di;
      dr /= den;
      di /= -den;
      dlr = cr*dr-ci*di;
      dli = cr*di+ci*dr;
      double temp = p*dlr-q*dli;
      q = p*dli+q*dlr;
      p = temp;
      if (Maths::abs(dlr-1.0)+Maths::abs(dli) <= Maths::absEpsilon) break;
    }
    if (i >= MAXIT) throw("cf2 failed in besselJY");
    double gam = (p-f)/q;
    rjmu = std::sqrt(w/((p-f)*gam+q));
    rjmu = Maths::sign(rjmu,rjl);
    rymu = rjmu*gam;
    double rymup = rymu*(p+q/gam);
    ry1 = xmu*xi*rymu-rymup;
  }

  fact = rjmu/rjl;
  Jo = rjl1*fact;
  Jpo = rjp1*fact;
  for (int i = 1; i <= nl; i++) {
    double rytemp = (xmu+i)*xi2*ry1-rymu;
    rymu = ry1;
    ry1 = rytemp;
  }
  Yo = rymu;
  Ypo = nu*xi*rymu-ry1;
  xJY = x;
  nuJY = nu;
}

void Bessel::besselIK(const double nu, const double x)
{ // Sets Io, Ko, Ipo, and Kpo respectively to the Bessel functions Inu(x), Knu(x) and their derivatives Inu'(x),
  // Knu'(x), for positive x and for xnu = \nu >= 0. The relative accuracy is within one or two significant digits of
  // Maths::absEpsilon. Fmin is a number close to the machine's smallest floating-point number.
  double rk1, rkmu, rkmup;
  int i;
  if (x <= 0.0 || nu < 0.0) throw("bad arguments in besselIK");
  // nl is the number of downward recurrences of the I's and upward recurrences of K's. xmu lies between -1/2 and 1/2.
  int nl = static_cast<int>(nu+0.5);
  double xmu = nu-nl;
  double xmu2 = xmu*xmu;
  double xi = 1.0/x;
  double xi2 = 2.0*xi;
  double h = nu*xi;
  if (h < Fmin) h = Fmin;
  double b = xi2*nu;
  double d = 0.0;
  double c = h;
  for (i = 0; i < MAXIT; i++) {
    b += xi2;
    d = 1.0/(b+d);
    c = b+1.0/c;
    double del = c*d;
    h = del*h;
    if (Maths::abs(del-1.0) <= Maths::absEpsilon) break;
  }
  if (i >= MAXIT) throw("x too large in besselIK; try asymptotic expansion");
  double ril = Fmin;
  double ripl = h*ril;
  double ril1 = ril;
  double rip1 = ripl;
  double fact = nu*xi;
  for (int l = nl-1; l >= 0; l--) {
    double ritemp = fact*ril+ripl;
    fact -= xi;
    ripl = fact*ritemp+ril;
    ril = ritemp;
  }
  double f = ripl/ril;

  if (x < XMIN) {
    double x2 = 0.5*x;
    double pimu = Maths::pi*xmu;
    fact = (Maths::abs(pimu) < Maths::absEpsilon ? 1.0 : pimu/std::sin(pimu));
    d = -std::log(x2);
    double e = xmu*d;
    double fact2 = (Maths::abs(e) < Maths::absEpsilon ? 1.0 : std::sinh(e)/e);
    double xx = 8.0*xmu*xmu-1.0;
    double gam1 = chebevshev(c1, NUSE1, xx);
    double gam2 = chebevshev(c2, NUSE2, xx);
    double ff = fact*(gam1*std::cosh(e)+gam2*fact2*d);
    double sum = ff;
    e = std::exp(e);
    double p = 0.5*e/(gam2-xmu*gam1);
    double q = 0.5/(e*(gam2+xmu*gam1));
    c = 1.0;
    d = x2*x2;
    double sum1 = p;
    for (i = 1; i <= MAXIT; i++) {
      ff = (i*ff+p+q)/(i*i-xmu2);
      c *= (d/i);
      p /= (i-xmu);
      q /= (i+xmu);
      double del = c*ff;
      sum += del;
      sum1 += c*(p-i*ff);
      if (Maths::abs(del) < Maths::abs(sum)*Maths::absEpsilon) break;
    }
    if (i > MAXIT) throw("bessk series failed to converge");
    rkmu = sum;
    rk1 = sum1*xi2;
  } else { // Evaluate CF2 by Steed's algorithm (Ch 5.2), which is OK because there can be no zero denominators.
    b = 2.0*(1.0+x);
    d = 1.0/b;
    double delh = d;
    h = delh;
    double q1 = 0.0; // Initializations for recurrence (6.6.35).
    double q2 = 1.0;
    double a1 = 0.25-xmu2;
    double q = c = a1; // First term in equation (6.6.34).
    double a = -a1;
    double s = 1.0+q*delh;
    for (i = 1; i < MAXIT; i++) {
      a -= 2*i;
      c = -a*c/(i+1.0);
      double qnew = (q1-b*q2)/a;
      q1 = q2;
      q2 = qnew;
      q += c*qnew;
      b += 2.0;
      d = 1.0/(b+a*d);
      delh = (b*d-1.0)*delh;
      h += delh;
      double dels = q*delh;
      s += dels;
      // Need only test convergence of sum since CF2 itself converges more quickly.
      if (Maths::abs(dels/s) <= Maths::absEpsilon) break;
    }
    if (i >= MAXIT) throw("besselIK: failure to converge in cf2");
    h = a1*h;
    rkmu = std::sqrt(Maths::pi/(2.0*x))*std::exp(-x)/s;
    rk1 = rkmu*(xmu+x+0.5-h)*xi;
  }

  rkmup = xmu*xi*rkmu-rk1;
  double rimu = xi/(f*rkmu-rkmup);
  Io = (rimu*ril1)/ril;
  Ipo = (rimu*rip1)/ril;
  for (int i = 1; i <= nl; i++) {
    double rktemp = (xmu+i)*xi2*rk1+rkmu;
    rkmu = rk1;
    rk1 = rktemp;
  }
  Ko = rkmu;
  Kpo = nu*xi*rkmu-rk1;
  xIK = x;
  nuIK = nu;
}

void Bessel::airy(const double x)
{ // Sets Aio, Bio, Aipo, and Bipo to the Airy functions Ai(x), Bi(x) and their derivatives Ai0(x), Bi0(x).
  static const double ONOVRT = 0.577350269189626;
  double absx = Maths::abs(x);
  double rootx = std::sqrt(absx);
  double z = (2./3.)*absx*rootx;
  if (x > 0.0) {
    besselIK(1./3., z);
    Aio = rootx*ONOVRT*Ko/Maths::pi;
    Bio = rootx*(Ko/Maths::pi+2.0*ONOVRT*Io);
    besselIK(2./3., z);
    Aipo = -x*ONOVRT*Ko/Maths::pi;
    Bipo = x*(Ko/Maths::pi+2.0*ONOVRT*Io);
  } else if (x < 0.0) {
    besselJY(1./3., z);
    Aio = 0.5*rootx*(Jo-ONOVRT*Yo);
    Bio = -0.5*rootx*(Yo+ONOVRT*Jo);
    besselJY(2./3., z);
    Aipo = 0.5*absx*(ONOVRT*Yo+Jo);
    Bipo = 0.5*absx*(ONOVRT*Jo-Yo);
  } else {
    Aio = 0.355028053887817;
    Bio = Aio/ONOVRT;
    Aipo = -0.258819403792807;
    Bipo = -Aipo/ONOVRT;
  }
  xAB = x;
}

void Bessel::spherical(const int n, const double x)
{ // Sets jo, yo, jpo, and ypo, respectively, to the spherical Bessel functions jn(x), yn(x), and their
  // derivatives jn'(x), yn'(x) for integer n (which is saved as njy).
  if (n < 0 || x <= 0.0) throw("bad arguments in spherical");
  double order = n+0.5;
  besselJY(order, x);
  double factor = 1.253314137315500251/std::sqrt(x); // sqrt(pi/2x))
  jo = factor*Jo;

  yo = factor*Yo;
  jpo = factor*Jpo-jo/(2.*x);
  ypo = factor*Ypo-yo/(2.*x);
  njy = n;
  xjy = xJY;
}

const double Bessel::c1[7] = {-1.142022680371168e0, 6.5165112670737e-3, 3.087090173086e-4, -3.4706269649e-6, 6.9437664e-9, 3.67795e-11, -1.356e-13};
const double Bessel::c2[8] = {1.843740587300905e0, -7.68528408447867e-2, 1.2719271366546e-3, -4.9717367042e-6, -3.31261198e-8, 2.423096e-10, -1.702e-13, -1.49e-15};

#endif /* end of include guard: BESSEL_H */
