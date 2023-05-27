#include "../include/NutMatrix.h"
#include "../include/MeanObliquity.h"
#include "../include/NutAngles.h"
#include "../include/R_x.h"
#include "../include/R_z.h"
#include "../include/matrix.h"

void NutMatrix(double Mjd_TT, double **NutMat) {

  // Mean obliquity of the ecliptic
  double eps = MeanObliquity(Mjd_TT);

  // Nutation in longitude and obliquity
  double dpsi, deps;

  NutAngles(Mjd_TT, dpsi, deps);
  double **mat1 = new double *[3];
  double **mat2 = new double *[3];
  double **mat3 = new double *[3];

  for (int i = 0; i < 3; i++) {
    mat1[i] = new double[3];
    mat2[i] = new double[3];
    mat3[i] = new double[3];
  }
  // Transformation from mean to true equator and equinox
  R_x(-eps - deps, mat1);
  R_z(-dpsi, mat2);
  R_x(eps, mat3);
  mult3(mat1, 3, 3, mat2, 3, 3, mat3, 3, 3, NutMat);

  for (int i = 0; i < 3; i++) {
    delete[] mat1[i];
    delete[] mat2[i];
    delete[] mat3[i];
  }
  delete[] mat1;
  delete[] mat2;
  delete[] mat3;
}
