/**
 * Performs a time update.
 *
 * @param P Input/Output matrix
 * @param fp Rows of P
 * @param cp Columns of P
 * @param Phi Input matrix
 * @param fphi Rows of Phi
 * @param cphi Columns of Phi
 * @param Qdt Input matrix
 * @param fq Rows of Qdt
 * @param cq Columns of Qdt
 */
void TimeUpdate(double **P, int fp, int cp, double **Phi, int fphi, int cphi,
                double **Qdt = nullptr, int fq = 0, int cq = 0);
