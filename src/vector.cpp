#include <cmath>
#include "../include/vector.h"


double norm(double *v, int n) {
    double suma = 0;
    int i;
    if (n <= 0) {
        throw "Empty vector";
    }
    for (i = 0; i < n; i++) {
        suma += v[i] * v[i];
    }
    return sqrt(suma);
}

double dot(double *v1, int n1,double *v2, int n2) {
    double suma = 0;
    int i;
    if (n1 < 0 || n2 < 0 || n1 != n2) {
        throw "Invalid vectors";
    }
    for (i = 0; i < n1; i++) {
        suma += v1[i] * v2[i];
    }
    return suma;
}

void cross(double *v1, double *v2, double *v3) {
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}


