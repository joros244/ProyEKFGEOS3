// Global variables
extern double **eopdata;

/*
 * Loads the global variable eopdata as a 2d array of size [19716x13].
 * @param path Path to the eopdata file
 */
void loadEOP(const char *path);

/*
 * Deletes the global variable eopdata.
 */
void deleteEOP();
