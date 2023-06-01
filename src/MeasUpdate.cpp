#include "../include/MeasUpdate.h"
#include "../include/matrix.h"

void MeasUpdate(double **x, double z, double g, double s, double *G, double **P,
                int n, double **K) {

  double **Inv_W = new double *[1];
  Inv_W[0] = new double[1];
  Inv_W[0][0] = s * s; // Inverse weight (measurement covariance)

  //  Kalman gain
  double **GT = new double *[n];
  for (int i = 0; i < n; i++) {
    GT[i] = new double[1];
  }

  double **GG = new double *[1];
  GG[0] = new double[n];
  for (int i = 0; i < n; i++) {
    GG[0][i] = G[i];
  }

  transpose(GG, GT, 1, n);

  double **temp = new double *[1];
  temp[0] = new double[1];

  mult3(GG, 1, n, P, n, n, GT, n, 1, temp);

  double **inv_temp = new double *[1];
  inv_temp[0] = new double[1];
  temp[0][0] += Inv_W[0][0];
  inverMat(temp, 1, 1, inv_temp);

  mult3(P, n, n, GT, n, 1, inv_temp, 1, 1, K);

  // State update
  for (int i = 0; i < n; i++) {
    x[i][0] += K[i][0] * (z - g);
  }

  // Covariance update
  double **temp2 = new double *[n];
  for (int i = 0; i < n; i++) {
    temp2[i] = new double[n];
  }
  mult(K, n, 1, GG, 1, n, temp2);

  int k;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      k = (i == j) ? 1.0 : 0.0;
      temp2[i][j] = k - temp2[i][j];
    }
  }
  double **Pcp = new double *[n];
  for (int i = 0; i < n; i++) {
    Pcp[i] = new double[n];
  }
  mult(temp2, n, n, P, n, n, Pcp);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      P[i][j] = Pcp[i][j];
    }
  }

  for (int i = 0; i < n; i++) {
    delete[] GT[i];
    delete[] Pcp[i];
    delete[] temp2[i];
  }
  delete[] temp[0];
  delete[] GG[0];
  delete[] inv_temp[0];
  delete[] Inv_W[0];
  delete[] temp;
  delete[] inv_temp;
  delete[] Inv_W;
  delete[] GT;
  delete[] GG;
  delete[] temp2;
  delete[] Pcp;
}
