/**
 * Cimputes the azimuth, elevation, and partials from local tangent coordinates
 *
 *
 * @param s Topocentric local tangent coordinates (East-North-Zenith frame)
 * @param Az Azimuth angle in radians
 * @param Elevation in radians
 * @param dAds Partials of azimuth w.r.t. s
 * @param dEds Partials of elevation w.r.t. s
 *
 */
void AzElPa(double *s, double &Az, double &El, double *dAds, double *dEds);
