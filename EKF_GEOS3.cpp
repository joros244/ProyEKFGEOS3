//--------------------------------------------------------------------------
//
// Initial Orbit Determination using Gauss and Extended Kalman Filter methods
//
// References:
//   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
//   Applications", Springer Verlag, Heidelberg, 2000.
//
//   D. Vallado, "Fundamentals of Astrodynamics and Applications",
//   4th Edition, 2013.
//
//   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
//
// Last modified:   2020/03/16   Meysam Mahooti
//--------------------------------------------------------------------------
#include "include/Accel.h"
#include "include/AzElPa.h"
#include "include/DEInteg.h"
#include "include/IERS.h"
#include "include/LTC.h"
#include "include/MeasUpdate.h"
#include "include/Position.h"
#include "include/R_z.h"
#include "include/SAT_Const.h"
#include "include/TimeUpdate.h"
#include "include/VarEqn.h"
#include "include/anglesdr.h"
#include "include/global.h"
#include "include/gmst.h"
#include "include/matrix.h"
#include "include/mjday.h"
#include "include/timediff.h"
#include "include/vector.h"
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>
using namespace std;

int main() {
  string path = "data/eop19620101.txt";
  loadEOP(path.c_str());
  string pathCS = "data/egm.txt";
  loadCS(pathCS.c_str());

  int nobs = 18;
  int Y, M, D, h, m;
  double s, az, el, Dist;
  double obs[18][4];

  // read observations
  FILE *fp;
  fp = fopen("data/GEOS3.txt", "r");
  for (int i = 0; i < nobs; i++) {

    fscanf(fp, "%i/%i/%i %i:%i:%lf %lf %lf %lf", &Y, &M, &D, &h, &m, &s, &az,
           &el, &Dist);

    obs[i][0] = mjday(Y, M, D, h, m, s);
    obs[i][1] = Rad * az;
    obs[i][2] = Rad * el;
    obs[i][3] = 1.0e3 * Dist;
  }
  fclose(fp);

  double sigma_range = 92.5;      // [m]
  double sigma_az = 0.0224 * Rad; // [rad]
  double sigma_el = 0.0139 * Rad; // [rad]

  // Kaena Point station
  double lat = Rad * 21.5748;     // [rad]
  double lon = Rad * (-158.2706); // [rad]
  double alt = 300.20;            // [m]

  double *Rs = new double[3];
  position(lon, lat, alt, Rs);

  double Mjd1 = obs[0][0];
  double Mjd2 = obs[8][0];
  double Mjd3 = obs[17][0];

  double *r2 = new double[3];
  double *v2 = new double[3];
  double *Rs1 = new double[3];
  double *Rs2 = new double[3];
  double *Rs3 = new double[3];

  for (int i = 0; i < 3; i++) {
    Rs1[i] = Rs[i];
    Rs2[i] = Rs[i];
    Rs3[i] = Rs[i];
  }

  anglesdr(eopdata, obs[0][1], obs[8][1], obs[17][1], obs[0][2], obs[8][2],
           obs[17][2], Mjd1, Mjd2, Mjd3, Rs1, Rs2, Rs3, r2, v2);

  double **Y0_apr = new double *[6];

  for (int i = 0; i < 6; i++) {
    Y0_apr[i] = new double[1];
    if (i <= 2) {
      Y0_apr[i][0] = r2[i];

    } else {
      Y0_apr[i][0] = v2[i - 3];
    }
  }

  double Mjd0 = mjday(1995, 1, 29, 02, 38, 00.0);

  double Mjd_UTC = obs[8][0];
  double UT1_UTC, TAI_UTC, x_pole, y_pole;
  IERS(eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
  double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
  timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

  AuxParam.Mjd_TT = Mjd_UTC + TT_UTC / 86400;
  AuxParam.n = 10;
  AuxParam.m = 10;

  n_eqn = 6;

  DEInteg(&Accel, 0.0, -(obs[8][0] - Mjd0) * 86400.0, 1.0e-13, 1.0e-6, 6,
          Y0_apr);

  double **P = new double *[6];
  double **Phi = new double *[6];
  for (int i = 0; i < 6; i++) {
    P[i] = new double[6];
    Phi[i] = new double[6];
  }

  for (int i = 0; i < 3; i++) {
    P[i][i] = 1.0e8;
  }
  for (int i = 3; i < 6; i++) {
    P[i][i] = 1.0e3;
  }
  double **LT = new double *[3];
  double **U = new double *[3];
  double **ss = new double *[3];
  double **r = new double *[3];
  double **aux1 = new double *[3];
  double *dAds = new double[3];
  double *dEds = new double[3];
  double *dDds = new double[3];
  double *aux2 = new double[3];
  double *aux3 = new double[3];
  for (int i = 0; i < 3; i++) {
    LT[i] = new double[3];
    U[i] = new double[3];
    ss[i] = new double[1];
    r[i] = new double[1];
    aux1[i] = new double[1];
  }

  LTC(lon, lat, LT);

  double **yPhi = new double *[42];
  for (int i = 0; i < 42; i++) {
    yPhi[i] = new double[1];
  }

  // Measurement loop
  double t = 0.0;
  double t_old = 0.0;
  double **Y_old = new double *[6];
  double *dAdY = new double[6];
  double *dEdY = new double[6];
  double *dDdY = new double[6];
  double **K = new double *[6];
  for (int i = 0; i < 6; i++) {
    Y_old[i] = new double[1];
    K[i] = new double[1];
  }

  double Azim = 0.0, Elev = 0.0, Mjd_TT = 0.0, Mjd_UT1 = 0.0, theta = 0.0;
  for (int ji = 0; ji < nobs; ji++) {
    // Previous step
    t_old = t;
    for (int ii = 0; ii < 6; ii++) {
      Y_old[ii][0] = Y0_apr[ii][0];
    }
    // Time increment and propagation
    Mjd_UTC = obs[ji][0];           // Modified Julian Date
    t = (Mjd_UTC - Mjd0) * 86400.0; // Time since epoch [s]

    IERS(eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC / 86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    for (int ii = 0; ii < 6; ii++) {
      yPhi[ii][0] = Y_old[ii][0];
      for (int j = 0; j < 6; j++) {
        if (ii == j) {
          yPhi[6 * (j + 1) + ii][0] = 1.0;
        } else {
          yPhi[6 * (j + 1) + ii][0] = 0.0;
        }
      }
    }

    DEInteg(&VarEqn, 0.0, t - t_old, 1.0e-13, 1.0e-6, 42, yPhi);

    // Extract state transition matrices
    for (int j = 0; j < 6; j++) {
      for (int ii = 0; ii < 6; ii++) {
        Phi[ii][j] = yPhi[6 * (j + 1) + ii][0];
      }
    }

    DEInteg(&Accel, 0.0, t - t_old, 1.0e-13, 1.0e-6, 6, Y0_apr);

    // Topocentric coordinates
    theta = gmst(Mjd_UT1); // Earth rotation
    R_z(theta, U);
    r[0][0] = Y0_apr[0][0];
    r[1][0] = Y0_apr[1][0];
    r[2][0] = Y0_apr[2][0];
    mult(U, 3, 3, r, 3, 1, aux1);
    aux1[0][0] -= Rs[0];
    aux1[1][0] -= Rs[1];
    aux1[2][0] -= Rs[2];
    mult(LT, 3, 3, aux1, 3, 1, ss); // Topocentric position [m]
    transpose(ss, &aux2, 3, 1);

    // Time update
    TimeUpdate(P, 6, 6, Phi, 6, 6);

    // Azimuth and partials
    AzElPa(aux2, Azim, Elev, dAds, dEds); // Azimuth, Elevation
    mult3(&dAds, 1, 3, LT, 3, 3, U, 3, 3, &aux3);
    dAdY[0] = aux3[0];
    dAdY[1] = aux3[1];
    dAdY[2] = aux3[2];

    // Measurement update
    MeasUpdate(Y0_apr, obs[ji][1], Azim, sigma_az, dAdY, P, 6, K);

    // Elevation and partials
    r[0][0] = Y0_apr[0][0];
    r[1][0] = Y0_apr[1][0];
    r[2][0] = Y0_apr[2][0];
    mult(U, 3, 3, r, 3, 1, aux1);
    aux1[0][0] -= Rs[0];
    aux1[1][0] -= Rs[1];
    aux1[2][0] -= Rs[2];
    mult(LT, 3, 3, aux1, 3, 1, ss); // Topocentric position [m]
    transpose(ss, &aux2, 3, 1);
    AzElPa(aux2, Azim, Elev, dAds, dEds); // Azimuth, Elevation
    mult3(&dEds, 1, 3, LT, 3, 3, U, 3, 3, &aux3);
    dEdY[0] = aux3[0];
    dEdY[1] = aux3[1];
    dEdY[2] = aux3[2];

    // Measurement update
    MeasUpdate(Y0_apr, obs[ji][2], Elev, sigma_el, dEdY, P, 6, K);

    // Range and partials
    r[0][0] = Y0_apr[0][0];
    r[1][0] = Y0_apr[1][0];
    r[2][0] = Y0_apr[2][0];
    mult(U, 3, 3, r, 3, 1, aux1);
    aux1[0][0] -= Rs[0];
    aux1[1][0] -= Rs[1];
    aux1[2][0] -= Rs[2];
    mult(LT, 3, 3, aux1, 3, 1, ss); // Topocentric position [m]
    transpose(ss, &aux2, 3, 1);
    Dist = norm(aux2, 3);
    for (int ii = 0; ii < 3; ii++) {
      dDds[ii] = ss[ii][0] / Dist;
    }

    // Range
    mult3(&dDds, 1, 3, LT, 3, 3, U, 3, 3, &aux3);
    dDdY[0] = aux3[0];
    dDdY[1] = aux3[1];
    dDdY[2] = aux3[2];
    // Measurement update
    MeasUpdate(Y0_apr, obs[ji][3], Dist, sigma_range, dDdY, P, 6, K);
  }

  IERS(eopdata, obs[17][0], UT1_UTC, TAI_UTC, x_pole, y_pole);
  timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
  Mjd_TT = Mjd_UTC + TT_UTC / 86400;
  AuxParam.Mjd_UTC = Mjd_UTC;
  AuxParam.Mjd_TT = Mjd_TT;

  DEInteg(&Accel, 0.0, -(obs[17][0] - obs[0][0]) * 86400.0, 1.0e-13, 1.0e-6, 6,
          Y0_apr);

  double Y_true[6] = {5753.173e3, 2673.361e3,  3440.304e3,
                      4.324207e3, -1.924299e3, -5.728216e3};

  printf("\nError of Position Estimation\n");
  printf("dX %10.1f [m]\n", Y0_apr[0][0] - Y_true[0]);
  printf("dY %10.1f [m]\n", Y0_apr[1][0] - Y_true[1]);
  printf("dZ %10.1f [m]\n", Y0_apr[2][0] - Y_true[2]);
  printf("\nError of Velocity Estimation\n");
  printf("dVx %8.1f [m/s]\n", Y0_apr[3][0] - Y_true[3]);
  printf("dVy %8.1f [m/s]\n", Y0_apr[4][0] - Y_true[4]);
  printf("dVz %8.1f [m/s]\n", Y0_apr[5][0] - Y_true[5]);
  for (int i = 0; i < 6; i++) {
    delete[] Y0_apr[i];
    delete[] Y_old[i];
    delete[] P[i];
    delete[] Phi[i];
    delete[] K[i];
  }
  for (int i = 0; i < 42; i++) {
    delete[] yPhi[i];
  }
  for (int i = 0; i < 3; i++) {
    delete[] LT[i];
    delete[] U[i];
    delete[] ss[i];
    delete[] r[i];
    delete[] aux1[i];
  }
  deleteEOP();
  deleteCS();
  delete[] r2;
  delete[] v2;
  delete[] Rs;
  delete[] Rs1;
  delete[] Rs2;
  delete[] Rs3;
  delete[] Y0_apr;
  delete[] Y_old;
  delete[] P;
  delete[] U;
  delete[] ss;
  delete[] r;
  delete[] yPhi;
  delete[] Phi;
  delete[] LT;
  delete[] aux1;
  delete[] aux2;
  delete[] aux3;
  delete[] dAds;
  delete[] dEds;
  delete[] dDds;
  delete[] dDdY;
  delete[] dAdY;
  delete[] dEdY;
  delete[] K;
  return 0;
}
