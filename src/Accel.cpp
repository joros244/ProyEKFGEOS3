#include "../include/Accel.h"
#include "../include/AccelHarmonic.h"
#include "../include/GHAMatrix.h"
#include "../include/IERS.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/PrecMatrix.h"
#include "../include/SAT_Const.h"
#include "../include/global.h"
#include "../include/matrix.h"
#include "../include/timediff.h"
#include <iostream>
#include <string>
using namespace std;

void Accel(double x, double **Y, double **dY) {
  string path = "data/eop19620101.txt";
  loadEOP(path.c_str());

  double UT1_UTC, TAI_UTC, x_pole, y_pole;
  IERS(eopdata, AuxParam.Mjd_TT + x / 86400, UT1_UTC, TAI_UTC, x_pole, y_pole);
  double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
  timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

  double Mjd_UT1 = AuxParam.Mjd_TT + x / 86400 + (UT1_UTC - TT_UTC) / 86400.0;

  double **P = new double *[3];
  double **N = new double *[3];
  double **T = new double *[3];
  double **Pole = new double *[3];
  double **GHAMat = new double *[3];
  double **E = new double *[3];

  for (int i = 0; i < 3; i++) {
    P[i] = new double[3];
    N[i] = new double[3];
    T[i] = new double[3];
    Pole[i] = new double[3];
    GHAMat[i] = new double[3];
    E[i] = new double[3];
  }

  PrecMatrix(MJD_J2000, AuxParam.Mjd_TT + x / 86400, P);
  NutMatrix(AuxParam.Mjd_TT + x / 86400, N);
  mult(N, 3, 3, P, 3, 3, T);
  PoleMatrix(x_pole, y_pole, Pole);
  GHAMatrix(Mjd_UT1, GHAMat);
  mult3(Pole, 3, 3, GHAMat, 3, 3, T, 3, 3, E);

  double **a = new double *[3];
  double **r = new double *[3];

  for (int i = 0; i < 3; i++) {
    a[i] = new double[1];
    r[i] = new double[1];
    r[i][0] = Y[i][0];
  }
  // Acceleration due to harmonic gravity field
  AccelHarmonic(r, E, AuxParam.n, AuxParam.m, a);

  dY[0][0] = Y[3][0];
  dY[1][0] = Y[4][0];
  dY[2][0] = Y[5][0];
  dY[3][0] = a[0][0];
  dY[4][0] = a[1][0];
  dY[5][0] = a[2][0];

  for (int i = 0; i < 3; i++) {
    delete[] P[i];
    delete[] N[i];
    delete[] T[i];
    delete[] Pole[i];
    delete[] GHAMat[i];
    delete[] E[i];
    delete[] a[i];
    delete[] r[i];
  }
  delete[] P;
  delete[] N;
  delete[] T;
  delete[] Pole;
  delete[] GHAMat;
  delete[] E;
  delete[] a;
  delete[] r;
  deleteEOP();
}
