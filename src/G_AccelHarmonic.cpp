#include "../include/G_AccelHarmonic.h"
#include "../include/AccelHarmonic.h"
void G_AccelHarmonic(double **r, double **U, int n_max, int m_max, double **G) {

  double d = 1.0; // Position increment [m]

  double **dr = new double *[3];
  double **a1 = new double *[3];
  double **a2 = new double *[3];
  double **m1 = new double *[3];
  double **m2 = new double *[3];
  for (int i = 0; i < 3; i++) {
    dr[i] = new double[1];
    dr[i][0] = 0.0;
    a1[i] = new double[1];
    a2[i] = new double[1];
    m1[i] = new double[1];
    m2[i] = new double[1];
  }

  // Gradient
  for (int i = 0; i < 3; i++) {
    for (int w = 0; w < 3; w++) {
      dr[w][0] = 0.0;
    }
    // Set offset in i-th component of the position vector
    dr[i][0] = d;
    // Acceleration difference
    for (int k = 0; k < 3; k++) {
      m1[k][0] = r[k][0] + dr[k][0] / 2.0;
      m2[k][0] = r[k][0] - dr[k][0] / 2.0;
    }
    AccelHarmonic(m1, U, n_max, m_max, a1);
    AccelHarmonic(m2, U, n_max, m_max, a2);
    // Derivative with respect to i-th axis
    for (int j = 0; j < 3; j++) {
      G[j][i] = (a1[j][0] - a2[j][0]) / d;
    }
  }
  for (int i = 0; i < 3; i++) {
    delete[] dr[i];
    delete[] a1[i];
    delete[] a2[i];
    delete[] m1[i];
    delete[] m2[i];
  }
  delete[] dr;
  delete[] a1;
  delete[] a2;
  delete[] m1;
  delete[] m2;
}
