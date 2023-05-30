/**
 * Computes the gradient of the Earth's harmonic gravity field
 *
 * @param r 3x1 array representing the position of the satellite
 * @param U 3x3 matrix representing the rotation matrix
 * @param n_max The maximum degree
 * @param m_max The maximum order
 * @param G 3x3 matrix to store gradient
 */
void G_AccelHarmonic(double **r, double **U, int n_max, int m_max, double **G);
