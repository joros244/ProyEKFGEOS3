#include "../include/anglesdr.h"
#include "../include/GHAMatrix.h"
#include "../include/Geodetic.h"
#include "../include/IERS.h"
#include "../include/LTC.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/PrecMatrix.h"
#include "../include/SAT_Const.h"
#include "../include/doubler.h"
#include "../include/global.h"
#include "../include/matrix.h"
#include "../include/timediff.h"
#include "../include/vector.h"
#include <cmath>
#include <cstring>
#include <iostream>

using namespace std;

void anglesdr(double az1, double az2, double az3, double el1, double el2,
              double el3, double Mjd1, double Mjd2, double Mjd3, double *rsite1,
              double *rsite2, double *rsite3, double *r2, double *v2) {

  double magr1in = 1.1 * R_Earth;
  double magr2in = 1.11 * R_Earth;
  char direct = 'y';
  double tol = 1e-8 * R_Earth;
  double pctchg = 0.005;
  double t1 = (Mjd1 - Mjd2) * 86400.0;
  double t3 = (Mjd3 - Mjd2) * 86400.0;
  double los1[3] = {cos(el1) * sin(az1), cos(el1) * cos(az1), sin(el1)};
  double los2[3] = {cos(el2) * sin(az2), cos(el2) * cos(az2), sin(el2)};
  double los3[3] = {cos(el3) * sin(az3), cos(el3) * cos(az3), sin(el3)};

  double lon1, lat1, h1;
  Geodetic(rsite1, lon1, lat1, h1);
  double lon2, lat2, h2;
  Geodetic(rsite2, lon2, lat2, h2);
  double lon3, lat3, h3;
  Geodetic(rsite3, lon3, lat3, h3);

  double **mat1 = new double *[3];
  double **mat2 = new double *[3];
  double **mat3 = new double *[3];

  for (int i = 0; i < 3; i++) {
    mat1[i] = new double[3];
    mat2[i] = new double[3];
    mat3[i] = new double[3];
  }

  LTC(lon1, lat1, mat1);
  LTC(lon2, lat2, mat2);
  LTC(lon3, lat3, mat3);

  // body-fixed system

  double **mat1T = new double *[3];

  for (int i = 0; i < 3; i++) {
    mat1T[i] = new double[3];
  }
  transpose(mat1, mat1T, 3, 3);

  double **los11 = new double *[3];
  double **los22 = new double *[3];
  double **los33 = new double *[3];

  for (int i = 0; i < 3; i++) {
    los11[i] = new double;
    los11[i][0] = los1[i];
    los22[i] = new double;
    los22[i][0] = los2[i];
    los33[i] = new double;
    los33[i][0] = los3[i];
  }

  double **los11R = new double *[3];
  double **los22R = new double *[3];
  double **los33R = new double *[3];

  for (int i = 0; i < 3; i++) {
    los11R[i] = new double;
    los22R[i] = new double;
    los33R[i] = new double;
  }

  mult(mat1T, 3, 3, los11, 3, 1, los11R);
  mult(mat1T, 3, 3, los22, 3, 1, los22R);
  mult(mat1T, 3, 3, los33, 3, 1, los33R);

  //// mean of date system (J2000)
  double Mjd_UTC = Mjd1;
  string path = "data/eop19620101.txt";
  loadEOP(path.c_str());
  double UT1_UTC, TAI_UTC, x_pole, y_pole;
  IERS(eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
  double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
  timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
  double Mjd_TT = Mjd_UTC + TT_UTC / 86400;
  double Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

  double **P = new double *[3];
  double **N = new double *[3];
  double **T = new double *[3];

  for (int i = 0; i < 3; i++) {
    P[i] = new double[3];
    N[i] = new double[3];
    T[i] = new double[3];
  }
  PrecMatrix(MJD_J2000, Mjd_TT, P);
  NutMatrix(Mjd_TT, N);
  mult(N, 3, 3, P, 3, 3, T);

  double **Pole = new double *[3];
  double **Gha = new double *[3];
  double **E = new double *[3];

  for (int i = 0; i < 3; i++) {
    Pole[i] = new double[3];
    Gha[i] = new double[3];
    E[i] = new double[3];
  }
  PoleMatrix(x_pole, y_pole, Pole);
  GHAMatrix(Mjd_UT1, Gha);
  mult3(Pole, 3, 3, Gha, 3, 3, T, 3, 3, E);

  double **ET = new double *[3];
  double **rsite1M = new double *[3];
  double **rsite1R = new double *[3];

  for (int i = 0; i < 3; i++) {
    ET[i] = new double[3];
    rsite1M[i] = new double;
    rsite1M[i][0] = rsite1[i];
    rsite1R[i] = new double;
  }
  transpose(E, ET, 3, 3);
  mult(ET, 3, 3, los11R, 3, 1, los11);
  mult(ET, 3, 3, rsite1M, 3, 1, rsite1R);

  for (int i = 0; i < 3; i++) {
    rsite1[i] = rsite1R[i][0];
  }
  /////////////////
  Mjd_UTC = Mjd2;
  IERS(eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
  timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
  Mjd_TT = Mjd_UTC + TT_UTC / 86400;
  Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

  PrecMatrix(MJD_J2000, Mjd_TT, P);
  NutMatrix(Mjd_TT, N);
  mult(N, 3, 3, P, 3, 3, T);
  PoleMatrix(x_pole, y_pole, Pole);
  GHAMatrix(Mjd_UT1, Gha);
  mult3(Pole, 3, 3, Gha, 3, 3, T, 3, 3, E);

  transpose(E, ET, 3, 3);
  mult(ET, 3, 3, los22R, 3, 1, los22);

  double **rsite2M = new double *[3];
  double **rsite2R = new double *[3];

  for (int i = 0; i < 3; i++) {
    rsite2M[i] = new double;
    rsite2M[i][0] = rsite2[i];
    rsite2R[i] = new double;
  }
  mult(ET, 3, 3, rsite2M, 3, 1, rsite2R);

  for (int i = 0; i < 3; i++) {
    rsite2[i] = rsite2R[i][0];
  }
  ///////////////////
  Mjd_UTC = Mjd3;
  IERS(eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
  timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
  Mjd_TT = Mjd_UTC + TT_UTC / 86400;
  Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;

  PrecMatrix(MJD_J2000, Mjd_TT, P);
  NutMatrix(Mjd_TT, N);
  mult(N, 3, 3, P, 3, 3, T);
  PoleMatrix(x_pole, y_pole, Pole);
  GHAMatrix(Mjd_UT1, Gha);
  mult3(Pole, 3, 3, Gha, 3, 3, T, 3, 3, E);

  transpose(E, ET, 3, 3);
  mult(ET, 3, 3, los33R, 3, 1, los33);

  double **rsite3M = new double *[3];
  double **rsite3R = new double *[3];

  for (int i = 0; i < 3; i++) {
    rsite3M[i] = new double;
    rsite3M[i][0] = rsite3[i];
    rsite3R[i] = new double;
  }
  mult(ET, 3, 3, rsite3M, 3, 1, rsite3R);

  for (int i = 0; i < 3; i++) {
    rsite3[i] = rsite3R[i][0];
  }

  double magr1old = 99999999.9;
  double magr2old = 99999999.9;
  double magrsite1 = norm(rsite1, 3);
  double magrsite2 = norm(rsite2, 3);
  double magrsite3 = norm(rsite3, 3);

  double *los1V = new double[3];
  double *los2V = new double[3];
  double *los3V = new double[3];

  for (int i = 0; i < 3; i++) {
    los1V[i] = los11[i][0];
    los2V[i] = los22[i][0];
    los3V[i] = los33[i][0];
  }
  double cc1 = 2.0 * dot(los1V, 3, rsite1, 3);
  double cc2 = 2.0 * dot(los2V, 3, rsite2, 3);
  int ktr = 0;

  // for (int i = 0; i < 3; i++) {
  //   cout << "r1: " << rsite1[i] << endl;
  //   cout << "r2: " << rsite2[i] << endl;
  //   cout << "r3: " << rsite3[i] << endl;
  // }

  // for (int i = 0; i < 3; i++) {
  //   cout << "l1: " << los1V[i] << endl;
  //   cout << "l2: " << los2V[i] << endl;
  //   cout << "l3: " << los3V[i] << endl;
  // }

  double f, a, magr2, deltae32, g, magr1o, deltar1, q2, f1delr1, f2delr1,
      pf1pr1, pf2pr1, magr2o, deltar2, q3, f1delr2, f2delr2, pf1pr2, pf2pr2,
      delta, delta1, delta2, f1, f2, q1, magr1;
  double *r3 = new double[3];

  while (fabs(magr1in - magr1old) > tol || fabs(magr2in - magr2old) > tol) {
    ktr = ktr + 1;
    doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1V, los2V,
            los3V, rsite1, rsite2, rsite3, t1, t3, direct, r2, r3, f1, f2, q1,
            magr1, magr2, a, deltae32);

    f = 1.0 - a / magr2 * (1.0 - cos(deltae32));
    g = t3 - sqrt(pow(a, 3) / GM_Earth) * (deltae32 - sin(deltae32));
    for (int i = 0; i < 3; i++) {
      v2[i] = (r3[i] - f * r2[i]) / g;
    }

    magr1o = magr1in;
    magr1in = (1.0 + pctchg) * magr1in;
    deltar1 = pctchg * magr1in;
    doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1V, los2V,
            los3V, rsite1, rsite2, rsite3, t1, t3, direct, r2, r3, f1delr1,
            f2delr1, q2, magr1, magr2, a, deltae32);

    pf1pr1 = (f1delr1 - f1) / deltar1;
    pf2pr1 = (f2delr1 - f2) / deltar1;

    magr1in = magr1o;
    deltar1 = pctchg * magr1in;
    magr2o = magr2in;
    magr2in = (1.0 + pctchg) * magr2in;
    deltar2 = pctchg * magr2in;

    doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1V, los2V,
            los3V, rsite1, rsite2, rsite3, t1, t3, direct, r2, r3, f1delr2,
            f2delr2, q3, magr1, magr2, a, deltae32);

    pf1pr2 = (f1delr2 - f1) / deltar2;
    pf2pr2 = (f2delr2 - f2) / deltar2;

    magr2in = magr2o;
    deltar2 = pctchg * magr2in;

    delta = pf1pr1 * pf2pr2 - pf2pr1 * pf1pr2;
    delta1 = pf2pr2 * f1 - pf1pr2 * f2;
    delta2 = pf1pr1 * f2 - pf2pr1 * f1;

    deltar1 = -delta1 / delta;
    deltar2 = -delta2 / delta;

    magr1old = magr1in;
    magr2old = magr2in;

    magr1in = magr1in + deltar1;
    magr2in = magr2in + deltar2;
  }
  doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1V, los2V, los3V,
          rsite1, rsite2, rsite3, t1, t3, direct, r2, r3, f1, f2, q1, magr1,
          magr2, a, deltae32);

  f = 1.0 - a / magr2 * (1.0 - cos(deltae32));
  g = t3 - sqrt(pow(a, 3) / GM_Earth) * (deltae32 - sin(deltae32));
  for (int i = 0; i < 3; i++) {
    v2[i] = (r3[i] - f * r2[i]) / g;
  }

  for (int i = 0; i < 3; i++) {
    delete[] mat1[i];
    delete[] mat2[i];
    delete[] mat3[i];
    delete[] mat1T[i];
    delete[] los11[i];
    delete[] los22[i];
    delete[] los33[i];
    delete[] los11R[i];
    delete[] los22R[i];
    delete[] los33R[i];
    delete[] P[i];
    delete[] N[i];
    delete[] T[i];
    delete[] Pole[i];
    delete[] Gha[i];
    delete[] E[i];
    delete[] ET[i];
    delete[] rsite1M[i];
    delete[] rsite1R[i];
    delete[] rsite2M[i];
    delete[] rsite2R[i];
  }
  delete[] mat1;
  delete[] mat2;
  delete[] mat3;
  delete[] mat1T;
  delete[] los11;
  delete[] los22;
  delete[] los33;
  delete[] los11R;
  delete[] los22R;
  delete[] los33R;
  delete[] P;
  delete[] N;
  delete[] T;
  delete[] Pole;
  delete[] Gha;
  delete[] E;
  delete[] ET;
  delete[] rsite1M;
  delete[] rsite1R;
  delete[] rsite2M;
  delete[] rsite2R;
  delete[] los1V;
  delete[] los2V;
  delete[] los3V;
  delete[] r3;
  deleteEOP();
}
