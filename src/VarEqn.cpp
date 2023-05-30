#include "../include/GHAMatrix.h"
#include "../include/IERS.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/PrecMatrix.h"
#include "../include/SAT_Const.h"
#include "../include/global.h"
#include "../include/matrix.h"
#include "../include/timediff.h"
#include <string>
using namespace std;

void VarEqn(double x, double *yPhi, double **yPhip) {

  // global AuxParam eopdata

  double UT1_UTC, TAI_UTC, x_pole, y_pole;
  string path = "../include/eop19620101.txt";
  loadEOP(path.c_str());
  IERS(eopdata, AuxParam.Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
  double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
  timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
  double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

  // Transformation matrix
  double **P = new double *[3];
  double **N = new double *[3];
  double **T = new double *[3];
  double **E = new double *[3];
  double **Pole = new double *[3];
  double **GhaMat = new double *[3];
  for (int i = 0; i < 3; i++) {
    P[i] = new double[3];
    N[i] = new double[3];
    T[i] = new double[3];
    E[i] = new double[3];
    Pole[i] = new double[3];
    GhaMat[i] = new double[3];
  }

  PrecMatrix(MJD_J2000, AuxParam.Mjd_TT + x / 86400, P);
  NutMatrix(AuxParam.Mjd_TT + x / 86400, N);
  mult(N, 3, 3, P, 3, 3, T);
  PoleMatrix(x_pole, y_pole, Pole);
  GHAMatrix(Mjd_UT1, GhaMat);
  mult3(Pole, 3, 3, GhaMat, 3, 3, T, 3, 3, E);

  // State vector components
  double **Phi = new double *[6];
  for (int i = 0; i < 6; i++) {
    Phi[i] = new double[6];
  }

  double *r = new double[3];
  double *v = new double[3];
  for (int i = 0; i < 3; i++) {
    r[i] = yPhi[i];
    v[i] = yPhi[i + 3];
  }

  // State transition matrix
  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 6; i++) {
      Phi[i][j] = yPhi[6 * (j + 1) + i];
    }
  }

  // Acceleration and gradient
  //  a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
  // G = G_AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

  // Time derivative of state transition matrix
  double **dfdy = new double *[6];
  double **Phip = new double *[6];
  for (int i = 0; i < 6; i++) {
    dfdy[i] = new double[6];
    Phip[i] = new double[6];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dfdy[i][j] = 0.0; // dv/dr(i,j)
      // dfdy[i + 3][j] = G[i][j]; // da/dr(i,j)
      if (i == j) {
        dfdy[i][j + 3] = 1.0;
      } else {
        dfdy[i][j + 3] = 0.0; // dv/dv(i,j)
      }
      dfdy[i + 3][j + 3] = 0.0; // da/dv(i,j)
    }
  }

  mult(dfdy, 6, 6, Phi, 6, 6, Phip);

  // Derivative of combined state vector and state transition matrix
  for (int i = 0; i < 3; i++) {
    yPhip[i][0] = v[i]; // dr/dt(i)
    // yPhip[i + 3] = a[i]; // dv/dt(i)
  }

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      yPhip[6 * (j + 1) + i][0] = Phip[i][j]; // dPhi/dt(i,j) ATENCIÓN AQUÍ
    }
  }
  deleteEOP();
  for (int i = 0; i < 3; i++) {
    delete[] P[i];
    delete[] N[i];
    delete[] T[i];
    delete[] Pole[i];
    delete[] E[i];
    delete[] GhaMat[i];
  }
  for (int i = 0; i < 6; i++) {
    delete[] Phi[i];
    delete[] Phip[i];
    delete[] dfdy[i];
  }
  delete[] P;
  delete[] N;
  delete[] T;
  delete[] Pole;
  delete[] GhaMat;
  delete[] E;
  delete[] Phi;
  delete[] Phip;
  delete[] r;
  delete[] dfdy;
  delete[] v;
}
