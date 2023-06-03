/**
 * Computes the acceleration of an Earth orbiting satellite.
 *
 * @param x Modified Julian Date (TT)
 * @param Y Satellite state vector [6x1]
 * @param dY Result acceleration vector [6x1]
 *
 * @note Global variable eopdata must be loaded first
 *
 */
void Accel(double x, double **Y, double **dY);
