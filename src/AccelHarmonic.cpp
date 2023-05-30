#include "../include/AccelHarmonic.h"
#include "../include/Legendre.h"
#include "../include/SAT_Const.h"
#include "../include/global.h"
#include "../include/matrix.h"
#include "../include/vector.h"
#include <cmath>
#include <iostream>
#include <string>
using namespace std;

void AccelHarmonic(double **r, double **E, int n_max, int m_max, double **a) {

  string path = "data/egm.txt";
  loadCS(path.c_str());

  double gm = 398600.4415e9;  // [m^3/s^2]; JGM3/EGM96
  double r_ref = 6378.1363e3; // Radius Earth [m]; JGM3/EGM96

  // Body-fixed position
  double **r_bf = new double *[3];
  for (int i = 0; i < 3; i++) {
    r_bf[i] = new double[1];
  }
  mult(E, 3, 3, r, 3, 1, r_bf);

  double *r_bf_t = new double[3];
  for (int i = 0; i < 3; i++) {
    r_bf_t[i] = r_bf[i][0];
  }

  // Auxiliary quantities
  double d = norm(r_bf_t, 3); // distance
  double latgc = asin(r_bf[2][0] / d);
  double lon = atan2(r_bf[1][0], r_bf[0][0]);

  double **pnm = new double *[n_max + 1];
  double **dpnm = new double *[n_max + 1];
  for (int i = 0; i < n_max + 1; i++) {
    pnm[i] = new double[n_max + 1];
    dpnm[i] = new double[n_max + 1];
  }

  Legendre(n_max, m_max, latgc, pnm, dpnm);

  double dUdr = 0.0;
  double dUdlatgc = 0.0;
  double dUdlon = 0.0;
  double q3 = 0.0, q2 = q3, q1 = q2;
  double b1 = 0.0, b2 = 0.0, b3 = 0.0;
  for (int n = 0; n < n_max + 1; n++) {
    b1 = (-gm / pow(d, 2)) * pow((r_ref / d), n) * (n + 1);
    b2 = (gm / d) * pow((r_ref / d), n);
    b3 = (gm / d) * pow((r_ref / d), n);
    for (int m = 0; m <= n; m++) {
      q1 += pnm[n][m] * (Cnm[n][m] * cos(m * lon) + Snm[n][m] * sin(m * lon));
      q2 += dpnm[n][m] * (Cnm[n][m] * cos(m * lon) + Snm[n][m] * sin(m * lon));
      q3 +=
          m * pnm[n][m] * (Snm[n][m] * cos(m * lon) - Cnm[n][m] * sin(m * lon));
    }
    dUdr += q1 * b1;
    dUdlatgc += q2 * b2;
    dUdlon += q3 * b3;
    q3 = 0.0;
    q2 = q3;
    q1 = q2;
  }

  // Body-fixed acceleration
  double r2xy = pow(r_bf[0][0], 2) + pow(r_bf[1][0], 2);

  double ax =
      (1 / d * dUdr - r_bf[2][0] / (pow(d, 2) * sqrt(r2xy)) * dUdlatgc) *
          r_bf[0][0] -
      (1 / r2xy * dUdlon) * r_bf[1][0];
  double ay =
      (1 / d * dUdr - r_bf[2][0] / (pow(d, 2) * sqrt(r2xy)) * dUdlatgc) *
          r_bf[1][0] +
      (1 / r2xy * dUdlon) * r_bf[0][0];
  double az = 1 / d * dUdr * r_bf[2][0] + sqrt(r2xy) / pow(d, 2) * dUdlatgc;

  double **a_bf = new double *[3];
  for (int i = 0; i < 3; i++) {
    a_bf[i] = new double[1];
  }
  a_bf[0][0] = ax;
  a_bf[1][0] = ay;
  a_bf[2][0] = az;

  double **ET = new double *[3];
  for (int i = 0; i < 3; i++) {
    ET[i] = new double[3];
  }
  transpose(E, ET, 3, 3);
  // Inertial acceleration
  mult(ET, 3, 3, a_bf, 3, 1, a);

  for (int i = 0; i < 3; i++) {
    delete[] r_bf[i];
    delete[] a_bf[i];
    delete[] ET[i];
  }
  for (int i = 0; i < n_max + 1; i++) {
    delete[] pnm[i];
    delete[] dpnm[i];
  }
  delete[] pnm;
  delete[] a_bf;
  delete[] ET;
  delete[] dpnm;
  delete[] r_bf;
  delete[] r_bf_t;
  deleteCS();
}
