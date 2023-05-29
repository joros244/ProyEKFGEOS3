#include "../include/AzElPa.h"
#include "../include/SAT_Const.h"
#include "../include/vector.h"
#include <cmath>

void AzElPa(double *s, double &Az, double &El, double *dAds, double *dEds) {

  const double pi2 = 2.0 * pi;

  const double rho = sqrt(s[0] * s[0] + s[1] * s[1]);

  // Angles
  Az = atan2(s[0], s[1]);

  if (Az < 0.0) {
    Az = Az + pi2;
  }

  El = atan(s[2] / rho);

  // Partials
  const double a = dot(s, 3, s, 3);
  dAds[0] = s[1] / (rho * rho);
  dAds[1] = -s[0] / (rho * rho);
  dAds[2] = 0.0;
  dEds[0] = (-s[0] * s[2] / rho) / a;
  dEds[1] = (-s[1] * s[2] / rho) / a;
  dEds[2] = rho / a;
}
