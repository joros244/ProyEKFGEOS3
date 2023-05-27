/**
 * Solves the problem of orbit determination using three optical sightings.
 * 
 * @param az1    Azimuth at t1 [rad]
 * @param az2    Azimuth at t2 [rad]
 * @param az3    Azimuth at t3 [rad]
 * @param el1    Elevation at t1 [rad]
 * @param el2    Elevation at t2 [rad]
 * @param el3    Elevation at t3 [rad]
 * @param Mjd1   Modified Julian Date of t1
 * @param Mjd2   Modified Julian Date of t2
 * @param Mjd3   Modified Julian Date of t3
 * @param rsite1 ijk site1 position vector [m]
 * @param rsite2 ijk site2 position vector [m]
 * @param rsite3 ijk site3 position vector [m]
 * @param r      Pointer to the output ijk position vector at t2 [m]
 * @param v      Pointer to the output ijk velocity vector at t2 [m/s]
 */
void anglesdr(double az1, double az2, double az3, double el1, double el2, 
		double el3, double Mjd1, double Mjd2, double Mjd3, 
		double  rsite1, double rsite2, double rsite3, double *r, 
		double *v );
