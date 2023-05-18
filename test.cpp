#include "include/R_x.h"
#include "include/matrix.h"
#include "include/mjday.h"
#include <assert.h>
#include <cmath>
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
  double res4[3][3] = {};
  R_x(M_PI / 3, res4);
  double res5[3][3] = {
      {1, 0, 0}, {0, 0.5, sqrt(3) / 2}, {0, -sqrt(3) / 2, 0.5}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabs(res4[i][j] - res5[i][j]) < pow(10, -12));
    }
  }
  cout << "R_x test passed" << endl;
  // END R_X TEST

  // BEGIN MJDAY TEST
  assert(fabs(mjday(2023, 4, 27, 18, 17, 9.34) -
              60061.76191365718841552734375) < pow(10, -12));
  assert(fabs(mjday(2023, 4, 27) - 60061.0) < pow(10, -12));
  assert(fabs(mjday(0, 0, 0) + 678987.0) < pow(10, -12));
  cout << "Mjday test passed" << endl;
  // END MJDAY TEST
  cout << "All test passed" << endl;

  return 0;
}
