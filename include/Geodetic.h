/**
 * Geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 * from given position vector (r [m])
 *
 * @param r An array of size 3 representing the position vector.
 * @param lon Reference to a double  to store the longitude.
 * @param lat Reference to a double to store the latitude.
 * @param h Reference to a double to store the altitude.
 */
void Geodetic(double *r,double& lon,double& lat,double& h);
