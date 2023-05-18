#include "../include/matrix.h"
#include <iostream>
using namespace std;

void transpose(double **mat, double **res, int f, int c) {
  for (int i = 0; i < f; i++) {
    for (int j = 0; j < c; j++) {
      res[j][i] = mat[i][j];
    }
  }
}

void mult(double **mat1, int f1, int c1, double **mat2, int f2, int c2,
          double **res) {
  if (c1 != f2) {
    cout << "Incompatible matrices" << endl;
    return;
  }
  for (int i = 0; i < f1; i++) {
    for (int j = 0; j < c2; j++) {
      res[i][j] = 0;
      for (int k = 0; k < c1; k++) {
        res[i][j] += mat1[i][k] * mat2[k][j];
      }
    }
  }
}

void mult3(double **mat1, int f1, int c1, double **mat2, int f2, int c2,
          double **mat3, int f3, int c3,double **res) {
  if (c1 != f2 || c2!=f3) {
    cout << "Incompatible matrices" << endl;
    return;
  }
  
  double **resP = new double *[f1];
  for (int i = 0; i < f1; i++) {
    resP[i] = new double[c2];
  }
  mult(mat1, f1, c1, mat2, f2, c2, resP);
  mult(resP, f1, c2, mat3, f3, c3, res);

  for (int i = 0; i < f1; i++) {
    delete[] resP[i];
  }
  delete[] resP;
}
