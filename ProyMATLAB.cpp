#include <fstream>
#include <iostream>

using namespace std;

int main() {

  double Cnm[361][361] = {};
  double Snm[361][361] = {};

  FILE *fp = fopen("egm.txt", "r");
  for (int n = 0; n <= 360; n++) {
    for (int m = 0; m <= n; m++) {
      double temp[6];
      fscanf(fp, "%lf %lf %lf %lf %lf %lf", &temp[0], &temp[1], &temp[2],
             &temp[3], &temp[4], &temp[5]);
      Cnm[n][m] = temp[2];
      Snm[n][m] = temp[3];
    }
  }
  fclose(fp);

  FILE *fp2 = fopen("eop19620101.txt", "r");
  while (!feof(fp2)) {
    double fila[13] = {};
    fscanf(fp2, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &fila[0],
           &fila[1], &fila[2], &fila[3], &fila[4], &fila[5], &fila[6], &fila[7],
           &fila[8], &fila[9], &fila[10], &fila[11], &fila[12]);
    for (int i = 0; i < 13; i++) {
      cout << fila[i] << " ";
    }
    cout << endl;
  }
  fclose(fp2);

  int nobs = 0;
  int Y, M, D, h, m;
  double s, az, el, Dist;

  FILE *fp3;
  fp3 = fopen("GEOS3.txt", "r");
  while (nobs < 18 && !feof(fp3)) {

    fscanf(fp, "%i/%i/%i %i:%i:%lf %lf %lf %lf", &Y, &M, &D, &h, &m, &s, &az,
           &el, &Dist);

    cout << "Y: " << Y << endl;
    cout << "M: " << M << endl;
    cout << "D: " << D << endl;
    cout << "h: " << h << endl;
    cout << "m: " << m << endl;
    cout << "s: " << s << endl;
    cout << "az: " << az << endl;
    cout << "el: " << el << endl;
    cout << "Dist: " << Dist << endl;
    nobs++;
  }
  fclose(fp);

  return 0;
}
