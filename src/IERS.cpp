#include "../include/IERS.h"
#include "../include/matrix.h"
#include <math.h>

void IERS(double **eop, double Mjd_TT, double &UT1_UTC, double &TAI_UTC,
          double &x_pole, double &y_pole) {

  double Arcs = 3600.0 * 180.0 / M_PI; // Arcseconds per radian

  double mj = round(Mjd_TT);

  bool b = false;
  double *eopVec = new double[13];

  for (int i = 0; i < 19716; i++) {
    if (eop[i][3] == mj) {
      for (int j = 0; j < 13; j++) {
        eopVec[j] = eop[i][j];
      }
      b = true;
      break;
    }
  }

  // Setting of IERS Earth rotation parameters
  // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])

  if (b) {
    UT1_UTC = eopVec[6];       // UT1-UTC time difference [s]
    TAI_UTC = eopVec[12];      // TAI-UTC time difference [s]
    x_pole = eopVec[4] / Arcs; // Pole coordinate [rad]
    y_pole = eopVec[5] / Arcs; // Pole coordinate [rad]
  } else {                     // If !b => Linear indexing
    UT1_UTC = eop[0][6];       // UT1-UTC time difference [s]
    TAI_UTC = eop[0][12];      // TAI-UTC time difference [s]
    x_pole = eop[0][4] / Arcs; // Pole coordinate [rad]
    y_pole = eop[0][5] / Arcs; // Pole coordinate [rad]
  }

  delete[] eopVec;
  eopVec = nullptr;
}
