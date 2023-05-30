// Global variables
extern double **eopdata;
extern double **Cnm;
extern double **Snm;
struct Aux {
  double Mjd_TT = 0;
  double Mjd_UTC = 0;
  int n = 0;
  int m = 0;
};
extern Aux AuxParam;

/*
 * Loads the global variable eopdata as a 2d array of size [19716x13].
 * @param path Path to the eopdata file
 */
void loadEOP(const char *path);

/*
 * Deletes the global variable eopdata.
 */
void deleteEOP();
/*
 * Loads the global variables Cnm and Snm as a 2d array of size [361x361].
 * @param path Path to the egm file
 */
void loadCS(const char *path);

/*
 * Deletes the global variables Cnm,Snm.
 */
void deleteCS();
