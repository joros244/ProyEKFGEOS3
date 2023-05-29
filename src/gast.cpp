#include "../include/gast.h"
#include "../include/EqnEquinox.h"
#include "../include/gmst.h"
#include <cmath>

double gast(double Mjd_UT1) {
  return fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2 * M_PI);
}
