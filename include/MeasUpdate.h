/**
 * Performs the measurement update.
 *
 * @param x Input/Output vector
 * @param z Input value
 * @param g Input value
 * @param s Input value
 * @param G Input vector
 * @param P Input/Output matrix
 * @param n Size of P
 * @param K Output matrix
 *
 */
void MeasUpdate(double **x, double z, double g, double s, double *G, double **P,
                int n, double **K);
