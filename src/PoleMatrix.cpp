#include "../include/PoleMatrix.h"
#include "../include/R_x.h"
#include "../include/R_y.h"
#include "../include/matrix.h"

void PoleMatrix(double xp, double yp, double **PoleMat) {
  double **mat1 = new double *[3];
  double **mat2 = new double *[3];

  for (int i = 0; i < 3; i++) {
    mat1[i] = new double[3];
    mat2[i] = new double[3];
  }
  R_y(-xp, mat1);
  R_x(-yp, mat2);
  mult(mat1, 3, 3, mat2, 3, 3, PoleMat);

  for (int i = 0; i < 3; i++) {
    delete[] mat1[i];
    delete[] mat2[i];
  }
  delete[] mat1;
  delete[] mat2;
}
