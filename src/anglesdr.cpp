//---------------------------------------------------------------------------
//
//  anglesdr.m
//
//  this function solves the problem of orbit determination using three
//  optical sightings.
//
//  inputs:
//    az1      - azimuth at t1               rad
//    az2      - azimuth at t2               rad
//    az3      - azimuth at t3               rad
//    el1      - elevation at t1             rad
//    el2      - elevation at t2             rad
//    el3      - elevation at t3             rad
//    Mjd1     - Modified julian date of t1
//    Mjd2     - Modified julian date of t2
//    Mjd3     - Modified julian date of t3
//    rsite1   - ijk site1 position vector   m
//    rsite2   - ijk site2 position vector   m
//    rsite3   - ijk site3 position vector   m
//
//  outputs:
//    r        - ijk position vector at t2   m
//    v        - ijk velocity vector at t2   m/s
//
// Last modified:   2015/08/12   M. Mahooti
//
//---------------------------------------------------------------------------
#include "../include/anglesdr.h"
#include "../include/Geodetic.h"
#include "../include/IERS.h"
#include "../include/LTC.h"
#include "../include/PrecMatrix.h"
#include "../include/SAT_Const.h"
#include "../include/global.h"
#include "../include/matrix.h"
#include "../include/timediff.h"
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
  double **mat2T = new double *[3];
  double **mat3T = new double *[3];

  for (int i = 0; i < 3; i++) {
    mat1T[i] = new double[3];
    mat2T[i] = new double[3];
    mat3T[i] = new double[3];
  }
  transpose(mat1, mat1T, 3, 3);
  transpose(mat2, mat2T, 3, 3);
  transpose(mat3, mat3T, 3, 3);

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
  mult(mat2T, 3, 3, los22, 3, 1, los22R);
  mult(mat3T, 3, 3, los33, 3, 1, los33R);

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
    P[i] = new double;
    N[i] = new double;
    T[i] = new double;
  }
  PrecMatrix(MJD_J2000, Mjd_TT, P);
  // N = NutMatrix(Mjd_TT);
  // T = N * P;
  // E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;
  //
  // los1 = E'*los1;
  // rsite1 = E'*rsite1;
  //
  // Mjd_UTC = Mjd2;
  //[UT1_UTC, TAI_UTC, x_pole, y_pole] = IERS(eopdata, Mjd_UTC);
  //[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
  // Mjd_TT = Mjd_UTC + TT_UTC/86400;
  // Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
  //
  // P = PrecMatrix(MJD_J2000,Mjd_TT);
  // N = NutMatrix(Mjd_TT);
  // T = N * P;
  // E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;
  //
  // los2 = E'*los2;
  // rsite2 = E'*rsite2;
  //
  // Mjd_UTC = Mjd3;
  //[UT1_UTC, TAI_UTC, x_pole, y_pole] = IERS(eopdata, Mjd_UTC);
  //[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
  // Mjd_TT = Mjd_UTC + TT_UTC/86400;
  // Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
  //
  // P = PrecMatrix(MJD_J2000,Mjd_TT);
  // N = NutMatrix(Mjd_TT);
  // T = N * P;
  // E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;
  //
  // los3 = E'*los3;
  // rsite3 = E'*rsite3;
  //
  // magr1old  = 99999999.9;
  // magr2old  = 99999999.9;
  // magrsite1 = norm(rsite1);
  // magrsite2 = norm(rsite2);
  // magrsite3 = norm(rsite3);
  //
  // cc1 = 2.0*dot(los1,rsite1);
  // cc2 = 2.0*dot(los2,rsite2);
  // ktr = 0;
  //
  // while (abs(magr1in-magr1old) > tol | abs(magr2in-magr2old) > tol)
  //     ktr = ktr + 1;
  //     [r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler(
  //     cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
  //                     los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
  //
  //     f  = 1.0 - a/magr2*(1.0-cos(deltae32));
  //     g  = t3 - sqrt(a^3/GM_Earth)*(deltae32-sin(deltae32));
  //     v2 = (r3 - f*r2)/g;
  //
  //     magr1o = magr1in;
  //     magr1in = (1.0+pctchg)*magr1in;
  //     deltar1 = pctchg*magr1in;
  //     [r2,r3,f1delr1,f2delr1,q2,magr1,magr2,a,deltae32] = doubler(
  //     cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
  //                            los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
  //     pf1pr1 = (f1delr1-f1)/deltar1;
  //     pf2pr1 = (f2delr1-f2)/deltar1;
  //
  //     magr1in = magr1o;
  //     deltar1 = pctchg*magr1in;
  //     magr2o = magr2in;
  //     magr2in = (1.0+pctchg)*magr2in;
  //     deltar2 = pctchg*magr2in;
  //     [r2,r3,f1delr2,f2delr2,q3,magr1,magr2,a,deltae32] = doubler(
  //     cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
  //                            los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
  //     pf1pr2 = (f1delr2-f1)/deltar2;
  //     pf2pr2 = (f2delr2-f2)/deltar2;
  //
  //     magr2in = magr2o;
  //     deltar2 = pctchg*magr2in;
  //
  //     delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
  //     delta1 = pf2pr2*f1 - pf1pr2*f2;
  //     delta2 = pf1pr1*f2 - pf2pr1*f1;
  //
  //     deltar1 = -delta1/delta;
  //     deltar2 = -delta2/delta;
  //
  //     magr1old = magr1in;
  //     magr2old = magr2in;
  //
  //     magr1in = magr1in + deltar1;
  //     magr2in = magr2in + deltar2;
  //
  // end;
  //
  //[r2,r3,f1,f2,q1,magr1,magr2,a,deltae32] = doubler(
  // cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,...
  //												   los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct);
  //
  // f  = 1.0 - a/magr2*(1.0-cos(deltae32));
  // g  = t3 - sqrt(a^3/GM_Earth)*(deltae32-sin(deltae32));
  // v2 = (r3 - f*r2)/g;

  for (int i = 0; i < 3; i++) {
    delete[] mat1[i];
    delete[] mat2[i];
    delete[] mat3[i];
    delete[] mat1T[i];
    delete[] mat2T[i];
    delete[] mat3T[i];
    delete[] los11[i];
    delete[] los22[i];
    delete[] los33[i];
    delete[] los11R[i];
    delete[] los22R[i];
    delete[] los33R[i];
    delete[] P[i];
    delete[] N[i];
    delete[] T[i];
  }
  delete[] mat1;
  delete[] mat2;
  delete[] mat3;
  delete[] mat1T;
  delete[] mat2T;
  delete[] mat3T;
  delete[] los11;
  delete[] los22;
  delete[] los33;
  delete[] los11R;
  delete[] los22R;
  delete[] los33R;
  delete[] P;
  delete[] N;
  delete[] T;
  deleteEOP();
}
