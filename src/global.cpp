#include "../include/global.h"
#include <iostream>

using namespace std;

double **eopdata;
double **Cnm;
double **Snm;
int n_eqn;
Aux AuxParam;

void loadEOP(const char *path) {
  FILE *fp;
  fp = fopen(path, "r");

  if (fp == NULL) {
    cout << "Fail open eop19620101.txt file" << endl;
    exit(EXIT_FAILURE);
  }
  eopdata = new double *[19716];
  if (eopdata == NULL) {
    cout << "eopdata: memory not allocated" << endl;
    exit(EXIT_FAILURE);
  }

  for (int k = 0; k < 19716; k++) {
    eopdata[k] = new double[13];
    if (eopdata[k] == NULL) {
      cout << "eopdata[i]: memory not allocated" << endl;
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < 19716; i++) {
    fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &eopdata[i][0],
           &eopdata[i][1], &eopdata[i][2], &eopdata[i][3], &eopdata[i][4],
           &eopdata[i][5], &eopdata[i][6], &eopdata[i][7], &eopdata[i][8],
           &eopdata[i][9], &eopdata[i][10], &eopdata[i][11], &eopdata[i][12]);
  }
  fclose(fp);
}

void deleteEOP() {
  for (int i = 0; i < 19716; i++) {
    delete[] eopdata[i];
  }
  delete[] eopdata;
}

void loadCS(const char *path) {
  FILE *fp;
  fp = fopen(path, "r");

  if (fp == NULL) {
    cout << "Fail open egm.txt file" << endl;
    exit(EXIT_FAILURE);
  }
  Cnm = new double *[361];
  Snm = new double *[361];

  for (int k = 0; k < 361; k++) {
    Cnm[k] = new double[361];
    Snm[k] = new double[361];
  }

  int f, c;
  double aux1, aux2;

  for (int n = 0; n <= 360; n++) {
    for (int m = 0; m <= n; m++) {
      fscanf(fp, "%d%d%lf%lf%lf%lf", &f, &c, &Cnm[n][m], &Snm[n][m], &aux1,
             &aux2);
    }
  }
  fclose(fp);
}
void deleteCS() {
  for (int i = 0; i < 361; i++) {
    delete[] Cnm[i];
    delete[] Snm[i];
  }
  delete[] Cnm;
  delete[] Snm;
}
