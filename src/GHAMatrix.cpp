#include "../include/GHAMatrix.h"
#include "../include/R_z.h"
#include "../include/gast.h"

void GHAMatrix(double Mjd_UT1, double **GHAmat) { R_z(gast(Mjd_UT1), GHAmat); }
