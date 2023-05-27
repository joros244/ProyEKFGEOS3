#include <cmath>
#include "../include/R_y.h"

void R_y(double angle, double **res){

double C = cos(angle);
double S = sin(angle);

res[0][0] =   C;  res[0][1] = 0.0;  res[0][2] = -1.0*S;
res[1][0] = 0.0;  res[1][1] = 1.0;  res[1][2] =    0.0;
res[2][0] =   S;  res[2][1] = 0.0;  res[2][2] =      C;
}
