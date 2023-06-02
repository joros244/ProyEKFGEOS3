/**
 * Management of IERS time and polar motion data
 *
 * @param eop      Pointer to the eopdata
 * @param Mjd_TT   Modified julian date representing the input time
 * @param UT1_UTC  Output parameter to store the UT1-UTC time
 * @param TAI_UTC  Output parameter to store the TAI-UTC time
 * @param x_pole   Output parameter to store the x-component pole
 * @param y_pole   Output parameter to store the y-component pole
 * @param eop Eopdata matrix
 */
void IERS(double **eop, double Mjd_TT, double &UT1_UTC, double &TAI_UTC,
          double &x_pole, double &y_pole);
