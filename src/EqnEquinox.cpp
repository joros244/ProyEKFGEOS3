#include "../include/EqnEquinox.h"
#include "../include/MeanObliquity.h"
#include "../include/NutAngles.h"
#include <cmath>

double EqnEquinox(double Mjd_TT) {

  // Nutation in longitude and obliquity
  double dpsi, deps;
  NutAngles(Mjd_TT, dpsi, deps);

  // Equation of the equinoxes
  return dpsi * cos(MeanObliquity(Mjd_TT));
}
