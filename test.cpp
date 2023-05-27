#include "include/Geodetic.h"
#include "include/IERS.h"
#include "include/LTC.h"
#include "include/MeanObliquity.h"
#include "include/NutAngles.h"
#include "include/NutMatrix.h"
#include "include/PoleMatrix.h"
#include "include/Position.h"
#include "include/PrecMatrix.h"
#include "include/R_x.h"
#include "include/R_y.h"
#include "include/R_z.h"
#include "include/global.h"
#include "include/matrix.h"
#include "include/mjday.h"
#include "include/timediff.h"
#include "include/vector.h"
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <math.h>
using namespace std;

int main() {
  // BEGIN MATRIX TEST
  double **matrix = new double *[3];
  for (int i = 0; i < 3; i++) {
    matrix[i] = new double[3];
    for (int j = 0; j < 3; j++) {
      matrix[i][j] = i + 1;
    }
  }
  double **res = new double *[3];
  for (int i = 0; i < 3; i++) {
    res[i] = new double[3];
  }

  transpose((double **)matrix, (double **)res, 3, 3);

  // Test traspose
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(matrix[i][j] - res[j][i]) < pow(10, -12));
    }
  }

  double **res2 = new double *[3];
  for (int i = 0; i < 3; i++) {
    res2[i] = new double[3];
  }

  mult(matrix, 3, 3, res, 3, 3, res2);

  int resMult[3][3] = {{3, 6, 9}, {6, 12, 18}, {9, 18, 27}};

  // Test matrix multipliction
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(res2[i][j] - resMult[i][j]) < pow(10, -12));
    }
  }

  double **vec = new double *[3];
  for (int i = 0; i < 3; i++) {
    vec[i] = new double;
    vec[i][0] = i + 1;
  }

  double **vecRes = new double *[3];
  for (int i = 0; i < 3; i++) {
    vecRes[i] = new double;
  }

  mult(matrix, 3, 3, vec, 3, 1, vecRes);

  int resVec[3][1] = {{6}, {12}, {18}};

  // Test matrix*vector
  for (int i = 0; i < 3; i++) {
    assert(fabs(vecRes[i][0] - resVec[i][0]) < pow(10, -12));
  }

  double **res3 = new double *[3];
  for (int i = 0; i < 3; i++) {
    res3[i] = new double[3];
  }

  mult3(matrix, 3, 3, res, 3, 3, res2, 3, 3, res3);

  int resMult3[3][3] = {{126, 252, 378}, {252, 504, 756}, {378, 756, 1134}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(res3[i][j] - resMult3[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] matrix[i];
    delete[] res[i];
    delete[] res2[i];
    delete[] vec[i];
    delete[] vecRes[i];
    delete[] res3[i];
  }
  delete[] matrix;
  delete[] res;
  delete[] res2;
  delete[] vec;
  delete[] vecRes;
  delete[] res3;
  cout << "Matrix test passed" << endl;

  // END MATRIX TEST

  // BEGIN R_X TEST
  double **res4 = new double *[3];
  for (int i = 0; i < 3; i++) {
    res4[i] = new double[3];
  }
  R_x(M_PI / 3, res4);
  double res5[3][3] = {
      {1, 0, 0}, {0, 0.5, sqrt(3) / 2}, {0, -sqrt(3) / 2, 0.5}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(res4[i][j] - res5[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] res4[i];
  }
  delete[] res4;
  cout << "R_x test passed" << endl;
  // END R_X TEST

  // BEGIN VECTOR TEST
  double *vecTest1 = new double[3];
  vecTest1[0] = 1;
  vecTest1[1] = 2;
  vecTest1[2] = 3;
  double *vecTest2 = new double[3];
  vecTest2[0] = 0;
  vecTest2[1] = 1;
  vecTest2[2] = 0;
  assert(fabs(norm(vecTest1, 3) - sqrt(14.0)) < pow(10, -12));
  assert(fabs(dot(vecTest1, 3, vecTest2, 3) - 2.0) < pow(10, -12));
  double *vecTest3 = (double *)malloc(3);
  cross(vecTest1, vecTest2, vecTest3);
  assert(fabs(vecTest3[0] + 3.0) < pow(10, -12));
  assert(fabs(vecTest3[1] - 0.0) < pow(10, -12));
  assert(fabs(vecTest3[2] - 1.0) < pow(10, -12));
  cout << "Vector test passed" << endl;
  delete[] vecTest1;
  delete[] vecTest2;
  delete[] vecTest3;
  // END VECTOR TEST

  // BEGIN MJDAY TEST
  assert(fabs(mjday(2023, 4, 27, 18, 17, 9.34) -
              60061.76191365718841552734375) < pow(10, -12));
  assert(fabs(mjday(2023, 4, 27) - 60061.0) < pow(10, -12));
  assert(fabs(mjday(0, 0, 0) + 678987.0) < pow(10, -12));
  cout << "Mjday test passed" << endl;
  // END MJDAY TEST

  // BEGIN POSITION TEST
  double *p = new double[3];
  p[0] = 0.0;
  p[1] = 0.0;
  p[2] = 0.0;
  position(0.0, 0.0, 0.0, p);
  double *resP = new double[3];
  resP[0] = 6378137;
  resP[1] = 0;
  resP[2] = 0;
  for (int i = 0; i < 3; i++) {
    assert(fabs(p[i] - resP[i]) < pow(10, -12));
  }
  position(1.0, 2.0, 3.0, p);
  resP[0] = -1438078.94344053;
  resP[1] = -2239675.25517784;
  resP[2] = 5776811.07900739;
  for (int i = 0; i < 3; i++) {
    assert(fabs(p[i] - resP[i]) < pow(10, -8));
  }
  delete[] p;
  delete[] resP;

  cout << "Position test passed" << endl;
  // END POSITION TEST

  // BEGIN GEODETIC TEST
  double *r = new double[3];
  r[0] = 1.0;
  r[1] = 2.0;
  r[2] = 3.0;
  double a = 0.0, b = 0.0, c = 0.0;
  Geodetic(r, a, b, c);
  assert(fabs(a - 1.10714871779409) < pow(10, -12));
  assert(fabs(b - 1.57074413624965) < pow(10, -12));
  assert(fabs(c + 6356749.31418683) < pow(10, -8));

  delete[] r;

  cout << "Geodetic test passed" << endl;
  // END GEODETIC TEST

  // BEGIN R_Y TEST
  double **ryTest = new double *[3];
  for (int i = 0; i < 3; i++) {
    ryTest[i] = new double[3];
  }
  R_y(M_PI / 3, ryTest);
  double resRyTest[3][3] = {
      {0.5, 0, -sqrt(3) / 2}, {0, 1, 0}, {sqrt(3) / 2, 0, 0.5}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(ryTest[i][j] - resRyTest[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] ryTest[i];
  }
  delete[] ryTest;
  cout << "R_y test passed" << endl;
  // END R_Y TEST

  // BEGIN R_Z TEST
  double **rzTest = new double *[3];
  for (int i = 0; i < 3; i++) {
    rzTest[i] = new double[3];
  }
  R_z(M_PI / 3, rzTest);
  double resRzTest[3][3] = {
      {0.5, sqrt(3) / 2, 0}, {-sqrt(3) / 2, 0.5, 0}, {0, 0, 1}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(rzTest[i][j] - resRzTest[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] rzTest[i];
  }
  delete[] rzTest;
  cout << "R_z test passed" << endl;
  // END R_Z TEST

  // BEGIN LTC TEST
  double **ltcTest = new double *[3];
  for (int i = 0; i < 3; i++) {
    ltcTest[i] = new double[3];
  }
  LTC(M_PI / 2, M_PI / 2, ltcTest);
  double resLtcTest[3][3] = {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(ltcTest[i][j] - resLtcTest[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] ltcTest[i];
  }
  delete[] ltcTest;
  cout << "LTC test passed" << endl;
  // END LTC TEST

  // BEGIN IERS TEST
  string path = "data/eop19620101.txt";
  loadEOP(path.c_str());
  double Mjtest = 31.45;
  double aa = 0.0, bb = 0.0, cc = 0.0, dd = 0.0;
  IERS(eopdata, Mjtest, aa, bb, cc, dd);
  assert(fabs(aa - 0.0326338) < pow(10, -12));
  assert(fabs(bb - 2.0) < pow(10, -12));
  assert(fabs(cc + 6.1571337500911e-8) < pow(10, -12));
  assert(fabs(dd - 1.03265314076331e-6) < pow(10, -12));
  deleteEOP();
  cout << "IERS test passed" << endl;
  // END IERS TEST

  // BEGIN TIMEDIFF TEST
  double u1 = 2.4, u2 = 1.2, u3, u4, u5, u6, u7;
  timediff(u1, u2, u3, u4, u5, u6, u7);
  assert(fabs(u3 - 1.2) < pow(10, -12));
  assert(fabs(u4 - 17.8) < pow(10, -12));
  assert(fabs(u5 - 20.2) < pow(10, -12));
  assert(fabs(u6 - 33.384) < pow(10, -12));
  assert(fabs(u7 + 17.8) < pow(10, -12));
  cout << "Timediff test passed" << endl;
  // END TIMEDIFF TEST

  // BEGIN PRECMATRIX TEST
  double **precTest = new double *[3];
  for (int i = 0; i < 3; i++) {
    precTest[i] = new double[3];
  }
  PrecMatrix(3.14, 2.78, precTest);
  double precResTest[3][3] = {
      {0.999999999999971, 0.000000220214640, 0.000000095832533},
      {-0.000000220214640, 0.999999999999976, -0.000000000000011},
      {-0.000000095832533, -0.000000000000011, 0.999999999999995}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(precTest[i][j] - precResTest[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] precTest[i];
  }
  delete[] precTest;
  cout << "PrecMatrix test passed" << endl;
  // END PRECMATRIX TEST

  // BEGIN MEANOBLIQUITY TEST
  assert(fabs(MeanObliquity(3.14) - 0.409413050674634) < pow(10, -12));
  cout << "MeanObliquity test passed" << endl;
  // END MEANOBLIQUITY TEST

  // BEGIN NUTANGLES TEST
  double dpsi, deps;
  NutAngles(3.14, dpsi, deps);
  assert(fabs(dpsi - 2.726049715404026e-05) < pow(10, -12));
  assert(fabs(deps - 3.873700829677442e-05) < pow(10, -12));
  cout << "NutAngles test passed" << endl;
  // END NUTANGLES TEST

  // BEGIN NUTMATRIX TEST
  double **nutMatTest = new double *[3];
  for (int i = 0; i < 3; i++) {
    nutMatTest[i] = new double[3];
  }
  NutMatrix(3.14, nutMatTest);
  double nutResTest[3][3] = {
      {0.999999999628433, -0.000025007543231, -0.000010851612159},
      {0.000025007122853, 0.999999998937039, -0.000038737143971},
      {0.000010852580868, 0.000038736872589, 0.999999999190838}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(nutMatTest[i][j] - nutResTest[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] nutMatTest[i];
  }
  delete[] nutMatTest;
  cout << "NutMatrix test passed" << endl;
  // END NUTMATRIX TEST

  // BEGIN POLEMATRIX TEST
  double **poleTest = new double *[3];
  for (int i = 0; i < 3; i++) {
    poleTest[i] = new double[3];
  }
  PoleMatrix(3.14, 2.78, poleTest);
  double poleResTest[3][3] = {
      {-0.999998731727540, 0.000563423816293, -0.001489663356476},
      {0, -0.935334586120739, -0.353764345301143},
      {-0.001592652916487, -0.353763896631566, 0.935333399861642}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(poleTest[i][j] - poleResTest[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] poleTest[i];
  }
  delete[] poleTest;
  cout << "PoleMatrix test passed" << endl;
  // END POLEMATRIX TEST
  cout << "All test passed" << endl;

  return 0;
}
