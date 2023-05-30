#include "../include/Legendre.h"
#include <cmath>
#include <iostream>
using namespace std;

void Legendre(int n, int m, double fi, double **pnm, double **dpnm) {

  pnm[0][0] = 1.0;
  dpnm[0][0] = 0.0;
  pnm[1][1] = sqrt(3.0) * cos(fi);
  dpnm[1][1] = -sqrt(3.0) * sin(fi);
  // diagonal coefficients

  for (int i = 1; i < n; i++) {
    pnm[i + 1][i + 1] =
        sqrt((2 * (i + 1) + 1.0) / (2.0 * (i + 1))) * cos(fi) * pnm[i][i];
  }
  for (int i = 1; i < n; i++) {
    dpnm[i + 1][i + 1] = sqrt((2 * (i + 1) + 1.0) / (2.0 * (i + 1))) *
                         ((cos(fi) * dpnm[i][i]) - (sin(fi) * pnm[i][i]));
  }
  // horizontal first step coefficients
  for (int i = 0; i < n; i++) {
    pnm[i + 1][i] = sqrt(2 * (i + 1) + 1.0) * sin(fi) * pnm[i][i];
  }
  for (int i = 0; i < n; i++) {
    dpnm[i + 1][i] = sqrt(2 * (i + 1) + 1.0) *
                     ((cos(fi) * pnm[i][i]) + (sin(fi) * dpnm[i][i]));
  }
  // horizontal second step coefficients
  int j = -1;
  int k = 1;
  while (1) {
    for (int i = k; i < n; i++) {
      pnm[i + 1][j + 1] =
          sqrt((2 * (i + 1) + 1.0) /
               (((i + 1) - (j + 1)) * ((i + 1) + (j + 1)))) *
          ((sqrt(2 * (i + 1) - 1.0) * sin(fi) * pnm[i][j + 1]) -
           (sqrt((((i + 1) + (j + 1) - 1.0) * ((i + 1) - (j + 1) - 1)) /
                 (2 * (i + 1) - 3.0)) *
            pnm[i - 1][j + 1]));
    }
    j = j + 1;
    k = k + 1;
    if (j > m) {
      break;
    }
  }
  j = -1;
  k = 1;
  while (1) {
    for (int i = k; i < n; i++) {
      dpnm[i + 1][j + 1] =
          sqrt((2 * (i + 1) + 1.0) /
               (((i + 1) - (j + 1)) * ((i + 1) + (j + 1)))) *
          ((sqrt(2 * (i + 1) - 1.0) * sin(fi) * dpnm[i][j + 1]) +
           (sqrt(2 * (i + 1) - 1.0) * cos(fi) * pnm[i][j + 1]) -
           (sqrt((((i + 1) + (j + 1) - 1.0) * ((i + 1) - (j + 1) - 1)) /
                 (2 * (i + 1) - 3)) *
            dpnm[i - 1][j + 1]));
    }
    j = j + 1;
    k = k + 1;
    if (j > m) {
      break;
    }
  }
}
