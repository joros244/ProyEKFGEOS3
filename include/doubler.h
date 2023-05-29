/**
 * Computes a few double precision values.
 *
 * @param cc1 Input value
 * @param cc2 Input value
 * @param magrsite1 Input value
 * @param magrsite2 Input value
 * @param magr1in Input value
 * @param magr2in Input value
 * @param los1 ijk vector
 * @param los2 ijk vector
 * @param los3 ijk vector
 * @param rsite1 ijk site1 position vector
 * @param rsite2 ijk site2 position vector
 * @param rsite3 ijk site3 position vector
 * @param t1 Time value
 * @param t3 Time value
 * @param direct Direction
 * @param r2 Position vector
 * @param r3 Position vector
 * @param f1 Output value
 * @param f2 Output value
 * @param q1 Output value
 * @param magr1 Output value
 * @param magr2 Output value
 * @param a Output value
 * @param deltae32 Output value
 */
void doubler(double cc1, double cc2, double magrsite1, double magrsite2,
             double magr1in, double magr2in, double *los1, double *los2,
             double *los3, double *rsite1, double *rsite2, double *rsite3,
             double t1, double t3, char direct, double *r2, double *r3,
             double &f1, double &f2, double &q1, double &magr1, double &magr2,
             double &a, double &deltae32);
