#include <cmath>
#include "../include/R_x.h"

void R_x(double angle, double **res){
  double C = cos(angle);
  double S = sin(angle);
  res[0][0] = 1.0; res[0][1] = 0.0; res[0][2] = 0.0;
  res[1][0] = 0.0; res[1][1] = C; res[1][2] = S;
  res[2][0] = 0.0; res[2][1] = -1.0*S; res[2][2] = C;
}
