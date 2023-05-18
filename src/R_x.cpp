#include <cmath>

void R_x(double a, double res[3][3]){
  double C = cos(a);
  double S = sin(a);
  res[0][0] = 1.0; res[0][1] = 0.0; res[0][2] = 0.0;
  res[1][0] = 0.0; res[1][1] = C; res[1][2] = S;
  res[2][0] = 0.0; res[2][1] = -1.0*S; res[2][2] = C;
}
