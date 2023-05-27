/**
 * Computes time differences
 *
 * @param UT1_UTC UT1-UTCtime difference [s]
 * @param TAI_UTC TAI-UTC time difference [s]
 * @param UT1_TAI UT1-TAI time difference [s]
 * @param UTC_GPS UTC-GPS time difference [s]
 * @param UT1_GPS UT1-GPS time difference [s]
 * @param TT_UTC  TT-UTC time difference [s]
 * @param GPS_UTC GPS-UTC time difference [s]
 */
void timediff(double UT1_UTC, double TAI_UTC, double &UT1_TAI, double &UTC_GPS,
              double &UT1_GPS, double &TT_UTC, double &GPS_UTC);
