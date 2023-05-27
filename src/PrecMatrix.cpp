#include "../include/PrecMatrix.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/SAT_Const.h"
#include "../include/matrix.h"

void PrecMatrix(double Mjd_1, double Mjd_2, double **PrecMat) {

  double T = (Mjd_1 - MJD_J2000) / 36525.0;
  double dT = (Mjd_2 - Mjd_1) / 36525.0;

  // Precession angles
  double zeta = ((2306.2181 + (1.39656 - 0.000139 * T) * T) +
                 ((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT) *
                dT / Arcs;

  double z = zeta + ((0.79280 + 0.000411 * T) + 0.000205 * dT) * dT * dT / Arcs;

  double theta = ((2004.3109 - (0.85330 + 0.000217 * T) * T) -
                  ((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT) *
                 dT / Arcs;

  // Precession matrix
  double **m1 = new double *[3];
  double **m2 = new double *[3];
  double **m3 = new double *[3];
  for (int i = 0; i < 3; i++) {
    m1[i] = new double[3];
    m2[i] = new double[3];
    m3[i] = new double[3];
  }
  R_z(-z, m1);
  R_y(theta, m2);
  R_z(-zeta, m3);
  mult3(m1, 3, 3, m2, 3, 3, m3, 3, 3, PrecMat);

  for (int i = 0; i < 3; i++) {
    delete[] m1[i];
    delete[] m2[i];
    delete[] m3[i];
  }
  delete[] m1;
  delete[] m2;
  delete[] m3;
}
