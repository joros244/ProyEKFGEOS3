/**
 * Transforms from Greenwich meridian system to local tangent coordinates.
 * 
 * @param lon  The longitude [rad]
 * @param lat  The latitude [rad]
 * @param M    Pointer to the result
 */
void LTC(double lon, double lat, double **M);
