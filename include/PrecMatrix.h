/**
 * Precession transformation of equatorial coordinates
 *
 * @param Epoch given (Modified Julian Date TT)
 * @param Epoch to precess to (Modified Julian Date TT)
 * @param PrecMat Precession transformation matrix
 */
void PrecMatrix(double Mjd_1, double Mjd_2, double **PrecMat);
