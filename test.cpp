#include "matriz.h"
#include <cmath>
#include <iostream>
using namespace std;

int main() {

  double **matriz = new double *[3];
  for (int i = 0; i < 3; i++) {
    matriz[i] = new double[3];
    for (int j = 0; j < 3; j++) {
      matriz[i][j] = i;
    }
  }
  double **res = new double *[3];
  for (int i = 0; i < 3; i++) {
    res[i] = new double[3];
  }
  transpuesta((double **)matriz, (double **)res, 3, 3);

  cout << "Entrada:" << endl;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << matriz[i][j] << " ";
    }
    cout << endl;
  }

  cout << "Transpuesta:" << endl;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << res[i][j] << " ";
    }
    cout << endl;
  }

  double **res2 = new double *[3];
  for (int i = 0; i < 3; i++) {
    res2[i] = new double[3];
  }

  mult(matriz, 3, 3, res, 3, 3, res2);

  cout << "Entrada*Entrada^t:" << endl;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << res2[i][j] << " ";
    }
    cout << endl;
  }

  cout << "Vector:" << endl;
  double **vec = new double *[3];
  for (int i = 0; i < 3; i++) {
    vec[i] = new double;
    vec[i][0] = i;
    cout << vec[i][0] << endl;
  }

  double **vecRes = new double *[3];
  for (int i = 0; i < 3; i++) {
    vecRes[i] = new double;
  }

  mult(matriz, 3, 3, vec, 3, 1, vecRes);

  cout << "Entrada*vector:" << endl;
  for (int i = 0; i < 3; i++) {
    cout << vecRes[i][0] << " ";
    cout << endl;
  }

  double **res3 = new double *[3];
  for (int i = 0; i < 3; i++) {
    res3[i] = new double[3];
  }

  mult3(matriz, 3, 3, res, 3, 3, res2, 3, 3, res3);

  cout << "(Entrada*Entrada^t)^2:" << endl;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << res3[i][j] << " ";
    }
    cout << endl;
  }

  double res4[3][3] = {};
  R_x(M_PI / 3, res4);

  cout << "Matriz R_x de pi/3:" << endl;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << res4[i][j] << " ";
    }
    cout << endl;
  }

  for (int i = 0; i < 3; i++) {
    delete[] matriz[i];
    delete[] res[i];
    delete[] res2[i];
    delete[] vec[i];
    delete[] vecRes[i];
    delete[] res3[i];
  }
  delete[] matriz;
  delete[] res;
  delete[] res2;
  delete[] vec;
  delete[] vecRes;
  delete[] res3;

  return 0;
}
