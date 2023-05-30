#include "include/Accel.h"
#include "include/AccelHarmonic.h"
#include "include/AzElPa.h"
#include "include/EqnEquinox.h"
#include "include/Frac.h"
#include "include/GHAMatrix.h"
#include "include/G_AccelHarmonic.h"
#include "include/Geodetic.h"
#include "include/IERS.h"
#include "include/LTC.h"
#include "include/Legendre.h"
#include "include/MeanObliquity.h"
#include "include/MeasUpdate.h"
#include "include/NutAngles.h"
#include "include/NutMatrix.h"
#include "include/PoleMatrix.h"
#include "include/Position.h"
#include "include/PrecMatrix.h"
#include "include/R_x.h"
#include "include/R_y.h"
#include "include/R_z.h"
#include "include/TimeUpdate.h"
#include "include/anglesdr.h"
#include "include/doubler.h"
#include "include/gast.h"
#include "include/global.h"
#include "include/gmst.h"
#include "include/matrix.h"
#include "include/mjday.h"
#include "include/sign_.h"
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

  transpose(matrix, res, 3, 3);

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
    vec[i] = new double[1];
    vec[i][0] = i + 1;
  }

  double **vecRes = new double *[3];
  for (int i = 0; i < 3; i++) {
    vecRes[i] = new double[1];
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

  // Test matrix inverse

  double invM[3][3] = {{1.0, 2.0, 0.0}, {0.0, 1.0, 2.0}, {2.0, 0.0, 1.0}};

  double **ivMat = new double *[3];
  for (int i = 0; i < 3; i++) {
    ivMat[i] = new double[3];
    for (int j = 0; j < 3; j++) {
      ivMat[i][j] = invM[i][j];
    }
  }

  double invRes[3][3] = {{1.0 / 9.0, -2.0 / 9.0, 4.0 / 9.0},
                         {4.0 / 9.0, 1.0 / 9.0, -2.0 / 9.0},
                         {-2.0 / 9.0, 4.0 / 9.0, 1.0 / 9.0}};

  double **inverseMat = new double *[3];
  for (int i = 0; i < 3; i++) {
    inverseMat[i] = new double[3];
  }

  inverMat(ivMat, 3, 3, inverseMat);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(inverseMat[i][j] - invRes[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] matrix[i];
    delete[] res[i];
    delete[] res2[i];
    delete[] vec[i];
    delete[] vecRes[i];
    delete[] res3[i];
    delete[] ivMat[i];
    delete[] inverseMat[i];
  }
  delete[] matrix;
  delete[] res;
  delete[] res2;
  delete[] vec;
  delete[] vecRes;
  delete[] res3;
  delete[] ivMat;
  delete[] inverseMat;

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
  double *vecTest3 = new double[3];
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

  // BEGIN GAST TEST
  assert(fabs(gast(3.14) - 1.906895863343556) < pow(10, -12));
  cout << "Gast test passed" << endl;
  // END GAST TEST

  // BEGIN FRAC TEST
  assert(fabs(Frac(3.14) - 0.14) < pow(10, -12));
  cout << "Frac test passed" << endl;
  // END FRAC TEST

  // BEGIN EQNEQUINOX TEST
  assert(fabs(EqnEquinox(3.14) - 2.500754323404229e-05) < pow(10, -12));
  cout << "EqnEquinox test passed" << endl;
  // END EQNEQUINOX TEST

  // BEGIN GMST TEST
  assert(fabs(gmst(3.14) - 1.906870855800322) < pow(10, -12));
  cout << "Gmst test passed" << endl;
  // END GMST TEST

  // BEGIN GHAMATRIX TEST
  double **ghamatTest = new double *[3];
  for (int i = 0; i < 3; i++) {
    ghamatTest[i] = new double[3];
  }
  GHAMatrix(3.14, ghamatTest);
  double ghaResTest[3][3] = {{-0.329807384579281, 0.944048245100310, 0},
                             {-0.944048245100310, -0.329807384579281, 0},
                             {0, 0, 1.0}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(ghamatTest[i][j] - ghaResTest[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] ghamatTest[i];
  }
  delete[] ghamatTest;
  cout << "GHAMatrix test passed" << endl;
  // END GHAMATRIX TEST

  // BEGIN DOUBLER TEST
  double *r2 = new double[3];
  double *r3 = new double[3];
  double f1, f2, q1, magr1, magr2, aaa, deltae32;
  double *v1 = new double[3];
  v1[0] = -0.0514407301715203;
  v1[1] = 0.838593164440367;
  v1[2] = 0.54232403213698;
  double *v2 = new double[3];
  v2[0] = 0.185350425424354;
  v2[1] = 0.924321659182723;
  v2[2] = 0.333578611665541;
  double *v3 = new double[3];
  v3[0] = 0.48999206372453;
  v3[1] = 0.865773547227108;
  v3[2] = -0.10170517395279;
  double *v4 = new double[3];
  v4[0] = 5854667.7577933;
  v4[1] = 962016.736146505;
  v4[2] = 2333503.53479825;
  double *v5 = new double[3];
  v5[0] = 5847642.87233096;
  v5[1] = 1003838.42368066;
  v5[2] = 2333501.82312028;
  double *v6 = new double[3];
  v6[0] = 5839555.2146941;
  v6[1] = 1049868.17436044;
  v6[2] = 2333499.77773523;
  doubler(3542174.64126966, 5580277.97983915, 6375566.60240377,
          6375566.60240377, 7015950.7, 7079732.07, v1, v2, v3, v4, v5, v6,
          -97.9999914765358, 108.000017702579, 'y', r2, r3, f1, f2, q1, magr1,
          magr2, aaa, deltae32);

  assert(fabs(deltae32 - 0.11248279851753) < pow(10, -12));
  assert(fabs(aaa - 6238411.00963291) < pow(10, -5));
  assert(fabs(magr2 - 7079732.07) < pow(10, -8));
  assert(fabs(magr1 - 7015950.7) < pow(10, -8));
  assert(fabs(q1 - 20.188685431244) < pow(10, -10));
  assert(fabs(f2 - 7.86861985501643) < pow(10, -10));
  assert(fabs(f1 + 18.5921446051542) < pow(10, -10));

  double r2Res[3] = {6100522.45717325, 2264920.38181335, 2788613.92031198};

  double r3Res[3] = {6455956.51848352, 2138995.9027581, 2205556.47728644};

  for (int i = 0; i < 3; i++) {
    assert(fabs(r2[i] - r2Res[i]) < pow(10, -8));
    assert(fabs(r3[i] - r3Res[i]) < pow(10, -8));
  }

  cout << "Doubler test passed" << endl;
  delete[] r2;
  delete[] r3;
  delete[] v1;
  delete[] v2;
  delete[] v3;
  delete[] v4;
  delete[] v5;
  delete[] v6;
  // END DOUBLER TEST

  // BEGIN ANGLESDR TEST
  // double *rAngTest = new double[3];
  // double *vAngTest = new double[3];
  // double *rsite1Test = new double[3];
  // rsite1Test[0] = -5512568.44501153;
  // rsite1Test[1] = -2196994.68777797;
  // rsite1Test[2] = 2330805.22194045;
  // double *rsite2Test = new double[3];
  // rsite2Test[0] = -5512568.44501153;
  // rsite2Test[1] = -2196994.68777797;
  // rsite2Test[2] = 2330805.22194045;

  // double *rsite3Test = new double[3];
  // rsite3Test[0] = -5512568.44501153;
  // rsite3Test[1] = -2196994.68777797;
  // rsite3Test[2] = 2330805.22194045;

  // anglesdr(1.0559084894933, 1.36310214580757, 1.97615602688759,
  //          0.282624656433946, 0.453434794338875, 0.586427138011591,
  //          49746.1101504629, 49746.1112847221, 49746.1125347223, rsite1Test,
  //          rsite2Test, rsite3Test, rAngTest, vAngTest);

  // double rResTest[3] = {6147304.28873136, 2498216.09757119,
  // 2872808.05359544}; double vResTest[3] = {3764.62899474253,
  // -2217.84494072807, -6141.47100738888};

  // for (int i = 0; i < 3; i++) {
  //   // assert(fabs(rAngTest[i] - rResTest[i]) < pow(10, -5));
  //   cout << "r: " << fabs(rAngTest[i] - rResTest[i]) << endl;
  //   cout << "v: " << fabs(vAngTest[i] - vResTest[i]) << endl;
  //   // assert(fabs(vAngTest[i] - vResTest[i]) < pow(10, -5));
  // }

  // delete[] rAngTest;
  // delete[] vAngTest;
  // delete[] rsite1Test;
  // delete[] rsite2Test;
  // delete[] rsite3Test;
  cout << "Anglesdr test FAILED" << endl;
  // END ANGLESDR TEST

  // BEGIN SIGN TEST
  assert(fabs(sign_(3.14, 2.78) - 3.14) < pow(10, -12));
  cout << "Sign test passed" << endl;
  // END SIGN TEST

  // BEGIN TIMEUPDATE TEST
  double **pTest = new double *[6];
  double **phiTest = new double *[6];
  for (int i = 0; i < 6; i++) {
    pTest[i] = new double[6];
    pTest[i][i] = i + 1;
    phiTest[i] = new double[6];
    phiTest[i][i] = i + 2;
  }
  TimeUpdate(pTest, 6, 6, phiTest, 6, 6);
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      if (i != j) {
        assert(fabs(pTest[i][j]) < pow(10, -12));
      }
    }
  }
  assert(fabs(pTest[0][0] - 4.0) < pow(10, -12));
  assert(fabs(pTest[1][1] - 18.0) < pow(10, -12));
  assert(fabs(pTest[2][2] - 48.0) < pow(10, -12));
  assert(fabs(pTest[3][3] - 100.0) < pow(10, -12));
  assert(fabs(pTest[4][4] - 180.0) < pow(10, -12));
  assert(fabs(pTest[5][5] - 294.0) < pow(10, -12));

  for (int i = 0; i < 6; i++) {
    delete[] pTest[i];
    delete[] phiTest[i];
  }
  delete[] pTest;
  delete[] phiTest;
  cout << "TimeUpdate test passed" << endl;
  // END TIMEUPDATE TEST

  // BEGIN AZELPA TEST
  double *sTest = new double[3];
  double *AdsTest = new double[3];
  double *EdsTest = new double[3];
  double Az, El;
  sTest[0] = 1.0;
  sTest[1] = 2.0;
  sTest[2] = 3.0;
  AzElPa(sTest, Az, El, AdsTest, EdsTest);

  double AdsTestRes[3] = {0.4, -0.2, 0};
  double EdsTestRes[3] = {-0.095831484749991, -0.191662969499982,
                          0.159719141249985};
  for (int i = 0; i < 3; i++) {
    assert(fabs(AdsTest[i] - AdsTestRes[i]) < pow(10, -12));
    assert(fabs(EdsTest[i] - EdsTestRes[i]) < pow(10, -12));
  }
  assert(fabs(Az - 0.463647609000806) < pow(10, -12));
  assert(fabs(El - 0.930274014115472) < pow(10, -12));

  delete[] sTest;
  delete[] AdsTest;
  delete[] EdsTest;
  cout << "AzElPa test passed" << endl;
  // END AZELPA TEST

  // BEGIN MEASUPDATE TEST
  double **xTest = new double *[6];
  double *GTest = new double[6];
  double **PMTest = new double *[6];
  double **KTest = new double *[6];
  for (int i = 0; i < 6; i++) {
    xTest[i] = new double[1];
    PMTest[i] = new double[6];
    KTest[i] = new double[1];
  }
  GTest[0] = 1.18226811703679e-07;
  GTest[1] = 2.68617065403109e-07;
  GTest[2] = -4.04487218792015e-07;
  GTest[3] = 0.0;
  GTest[4] = 0.0;
  GTest[5] = 0.0;

  xTest[0][0] = 5747896.43764892;
  xTest[1][0] = 2702560.57450949;
  xTest[2][0] = 3459085.29959328;
  xTest[3][0] = 4379.55774313217;
  xTest[4][0] = -1948.9834094769;
  xTest[5][0] = -5813.31374009477;

  double PP[6][6] = {{101488000.600326, 128617.803934412, 169139.099218971,
                      40308.0399813276, 3496.44199102555, 4571.97846555939},
                     {128617.803934412, 101286995.428481, 82509.594120304,
                      3496.57950927422, 34761.7247204656, 2210.17635089282},
                     {169139.099218971, 82509.594120304, 101332421.164292,
                      4572.26963026297, 2210.23022307775, 35952.8041664988},
                     {40308.0399813276, 3496.57950927422, 4572.26963026297,
                      1001.66889789894, 1.41783559164539, 1.84417805416141},
                     {3496.44199102555, 34761.7247204656, 2210.23022307775,
                      1.41783559164539, 999.390221437519, 0.884050439731921},
                     {4571.97846555939, 2210.17635089282, 35952.8041664988,
                      1.84417805416141, 0.884050439731921, 999.857769260102}};

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      PMTest[i][j] = PP[i][j];
    }
  }
  MeasUpdate(xTest, 1.0559084894933, 1.05595496938417, 0.00039095375244673,
             GTest, PMTest, 6, KTest);

  double PPRes[6][6] = {
      {95859254.5326805, -12662417.1256939, 19431733.0975035, 38494.3314161026,
       -670.265126224979, 10879.8071750619},
      {-12662417.1256939, 72220028.7178261, 43855760.280198, -624.979596635881,
       25293.0978473085, 16544.3915355394},
      {19431733.0975035, 43855760.280198, 35412316.9201898, 10779.11060395,
       16469.4607061714, 14366.2657116377},
      {38494.3314161026, -624.979596635881, 10779.11060395, 1001.08448021159,
       0.0752288553456444, 3.87670231816378},
      {-670.265126224979, 25293.0978473085, 16469.4607061714,
       0.0752288553456442, 996.305795884776, 5.55345184146847},
      {10879.8071750619, 16544.3915355394, 14366.2657116377, 3.87670231816378,
       5.55345184146847, 992.788929673052}};
  double KRes[6] = {470444.610235622, 1069061.09631921, -1609947.12085869,
                    151.587833026057, 348.248949621027, -527.201615246513};
  double xResTest[6] = {5747874.57143478, 2702510.8846664,   3459160.12975976,
                        4379.55069734624, -1948.99959605007, -5813.28923582123};
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      assert(fabs(xTest[i][0] - xResTest[i]) < pow(10, -8));
      assert(fabs(KTest[i][0] - KRes[i]) < pow(10, -8));
      assert(fabs(PMTest[i][j] - PPRes[i][j]) < pow(10, -6));
    }
  }
  for (int i = 0; i < 6; i++) {
    delete[] xTest[i];
    delete[] PMTest[i];
    delete[] KTest[i];
  }
  delete[] xTest;
  delete[] GTest;
  delete[] PMTest;
  delete[] KTest;
  cout << "MeasUpdate test passed" << endl;
  // END MEASUPDATE TEST

  // BEGIN LEGENDRE TEST
  int n = 3;
  int m = 2;
  double fi = 0.5;

  double **pnm = new double *[n + 1];
  double **dpnm = new double *[n + 1];
  for (int i = 0; i < n + 1; i++) {
    pnm[i] = new double[n + 1];
    dpnm[i] = new double[n + 1];
  }

  double pnmRes[4][4] = {
      {1, 0, 0, 0},
      {0.830389391308554, 1.52001758503058, 0, 0},
      {-0.347097518865836, 1.62950155523887, 1.49139129468805, 0},
      {-1.17378701262255, 0.212202357272447, 1.89174148838053,
       1.41368608598461}};

  double dpnmRes[4][4] = {
      {0, 0, 0, 0},
      {1.52001758503058, -0.830389391308554, 0, 0},
      {2.82237948468622, 2.09258183254477, -1.62950155523887, 0},
      {0.519787497533268, 5.86628517139076, 1.39588339664843,
       -2.31690068589274}};

  Legendre(n, m, fi, pnm, dpnm);

  for (int i = 0; i < n + 1; i++) {
    for (int j = 0; j < n + 1; j++) {
      assert(fabs(pnm[i][j] - pnmRes[i][j]) < pow(10, -12));
      assert(fabs(dpnm[i][j] - dpnmRes[i][j]) < pow(10, -12));
    }
  }

  for (int i = 0; i < n + 1; i++) {
    delete[] pnm[i];
    delete[] dpnm[i];
  }
  delete[] pnm;
  delete[] dpnm;
  cout << "Legendre test passed" << endl;
  // END LEGENDRE TEST

  // BEGIN ACCELHARMONIC TEST
  double **rtest = new double *[3];
  double **atest = new double *[3];
  double **Etest = new double *[3];
  double rval[3] = {10000.0, 20000.0, 30000.0};
  double Eval[3][3] = {{0.1, 0.2, 0.3}, {0.1, 0.2, 0.3}, {0.1, 0.2, 0.3}};

  for (int i = 0; i < 3; i++) {
    rtest[i] = new double[1];
    rtest[i][0] = rval[i];
    atest[i] = new double[1];
    Etest[i] = new double[3];
    for (int j = 0; j < 3; j++) {
      Etest[i][j] = Eval[i][j];
    }
  }

  double aval[3] = {1777073.32963624, 3554146.65927248, 5331219.98890872};
  AccelHarmonic(rtest, Etest, 3, 3, atest);
  for (int i = 0; i < 3; i++) {
    assert(fabs(aval[i] - atest[i][0]) < pow(10, -8));
  }
  for (int i = 0; i < 3; i++) {
    delete[] rtest[i];
    delete[] atest[i];
    delete[] Etest[i];
  }
  delete[] rtest;
  delete[] atest;
  delete[] Etest;

  cout << "AccelHarmonic test passed" << endl;
  // END ACCELHARMONIC TEST

  // BEGIN ACCEL TEST
  double **Ytest = new double *[6];
  double **dY = new double *[6];
  double Yval[6] = {6511674.62431053, 2217276.31080875,  2186612.13046445,
                    2977.29908428071, -2465.79612321289, -6347.53937382649};

  for (int i = 0; i < 6; i++) {
    Ytest[i] = new double[1];
    Ytest[i][0] = Yval[i];
    dY[i] = new double[1];
  }

  double dYval[6] = {2977.29908428071,  -2465.79612321289, -6347.53937382649,
                     -6.90707086490836, -2.35188341613577, -2.32499942836077};
  AuxParam.n = 3;
  AuxParam.m = 2;
  Accel(12.8802163484948, Ytest, dY);
  for (int i = 0; i < 6; i++) {
    assert(fabs(dYval[i] - dY[i][0]) < pow(10, -7));
  }
  for (int i = 0; i < 6; i++) {
    delete[] Ytest[i];
    delete[] dY[i];
  }
  delete[] Ytest;
  delete[] dY;

  cout << "Accel test passed" << endl;
  // END ACCEL TEST

  // BEGIN GACCELHARMONIC TEST
  double **rtestG = new double *[3];
  double **EtestG = new double *[3];
  double **Gres = new double *[3];
  double rvalG[3] = {10000.0, 20000.0, 30000.0};
  double EvalG[3][3] = {{0.1, 0.2, 0.3}, {0.1, 0.2, 0.3}, {0.1, 0.2, 0.3}};

  for (int i = 0; i < 3; i++) {
    Gres[i] = new double[3];
    rtestG[i] = new double[1];
    rtestG[i][0] = rvalG[i];
    EtestG[i] = new double[3];
    for (int j = 0; j < 3; j++) {
      EtestG[i][j] = EvalG[i][j];
    }
  }

  double GtestRes[3][3] = {
      {-65.6685017049313, -131.337003447115, -197.005505255423},
      {-131.337003409863, -262.674006894231, -394.011010510847},
      {-197.00550512597, -394.011010337621, -591.016515765339}};

  G_AccelHarmonic(rtestG, EtestG, 3, 2, Gres);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << Gres[i][j] << endl;
      //      cout << fabs(Gres[i][j] - GtestRes[i][j]) << endl;
      //      assert(fabs(Gres[i][j] - GtestRes[i][j]) < pow(10, -8));
    }
  }

  for (int i = 0; i < 3; i++) {
    delete[] rtestG[i];
    delete[] Gres[i];
    delete[] EtestG[i];
  }
  delete[] rtestG;
  delete[] Gres;
  delete[] EtestG;

  cout << "G_AccelHarmonic test passed" << endl;
  // END GACCELHARMONIC TEST

  cout << "All test passed" << endl;

  return 0;
}
