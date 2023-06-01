#include "../include/matrix.h"
#include <cmath>
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
           double **mat3, int f3, int c3, double **res) {
  if (c1 != f2 || c2 != f3) {
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

void inverMat(double **matrix, int rows, int columns, double **matrixInv) {

  double **matrixC = new double *[rows];
  double **matrixI = new double *[rows];
  for (int i = 0; i < rows; i++) {
    matrixC[i] = new double[columns];
    matrixI[i] = new double[columns];
  }

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      matrixC[i][j] = matrix[i][j];
      if (i == j) {
        matrixI[i][j] = 1;
      } else {
        matrixI[i][j] = 0;
      }
    }
  }

  int n = rows;

  for (int i = 0; i < n; i++) {
    int maxRow = i;
    for (int j = i + 1; j < n; j++) {
      if (fabs(matrixC[j][i]) > fabs(matrixC[maxRow][i])) {
        maxRow = j;
      }
    }

    if (maxRow != i) {
      swap(matrixC[i], matrixC[maxRow]);
      swap(matrixI[i], matrixI[maxRow]);
    }

    double pivot = matrixC[i][i];
    for (int j = 0; j < n; j++) {
      matrixC[i][j] = matrixC[i][j] / pivot;
      matrixI[i][j] = matrixI[i][j] / pivot;
    }

    for (int j = 0; j < n; j++) {
      if (j != i) {
        double factor = matrixC[j][i];
        for (int k = 0; k < n; k++) {
          matrixC[j][k] = matrixC[j][k] - factor * matrixC[i][k];
          matrixI[j][k] = matrixI[j][k] - factor * matrixI[i][k];
        }
      }
    }
  }

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      matrixInv[i][j] = matrixI[i][j];
    }
  }

  for (int i = 0; i < rows; i++) {
    delete[] matrixC[i];
    delete[] matrixI[i];
  }
  delete[] matrixC;
  delete[] matrixI;
}
