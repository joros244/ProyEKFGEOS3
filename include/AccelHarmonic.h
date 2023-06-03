/**
 * Computes the acceleration due to harmonic gravity field.
 *
 * @param r 3x1 array representing the position
 * @param E 3x3 matrix representing the transformation matrix
 * @param n_max The maximum degree
 * @param m_max The maximum order
 * @param a 3x1 array to store acceleration vector
 *
 * @note Global variables Cnm,Snm must be loaded first
 */
void AccelHarmonic(double **r, double **E, int n_max, int m_max, double **a);
