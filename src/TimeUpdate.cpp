#include "../include/TimeUpdate.h"
#include "../include/matrix.h"
#include <cstring>
#include <iostream>

using namespace std;

void TimeUpdate(double **P, int fp, int cp, double **Phi, int fphi, int cphi,
                double **Qdt, int fq, int cq) {

  double **aux1 = new double *[cphi];
  for (int i = 0; i < cphi; i++) {
    aux1[i] = new double[fphi];
  }
  transpose(Phi, aux1, fphi, cphi);

  double **aux2 = new double *[fphi];
  for (int i = 0; i < fphi; i++) {
    aux2[i] = new double[fphi];
  }
  mult3(Phi, fphi, cphi, P, fp, cp, aux1, cphi, fphi, aux2);

  for (int i = 0; i < fphi; i++) {
    for (int j = 0; j < fphi; j++) {
      if (!Qdt) {
        P[i][j] = aux2[i][j];
      } else {
        P[i][j] = aux2[i][j] + Qdt[i][j];
      }
    }
  }
  for (int i = 0; i < cphi; i++) {
    delete[] aux1[i];
  }
  for (int i = 0; i < fphi; i++) {
    delete[] aux2[i];
  }
  delete[] aux1;
  delete[] aux2;
}
