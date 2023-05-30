/**
 * Calculates the Legendre functions.
 *
 * @param n The maximum degree of the Legendre functions.
 * @param m The maximum order of the Legendre functions.
 * @param fi The latitude angle in radians.
 * @param pnm Result matrix [n+1 x n+1].
 * @param dpnm Result matrix [n+1 x n+1].
 *
 */
void Legendre(int n, int m, double fi, double **pnm, double **dpnm);
