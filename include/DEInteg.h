#include <functional>
using namespace std;
/**
 * Numerical integration methods for ODEs.
 *
 * @param func System of ODEs to be integrated
 * @param t       The initial time
 * @param tout    The final time
 * @param relerr  The relative error tolerance
 * @param abserr  The absolute error tolerance
 * @param n_eqn   The number of equations
 * @param y   Initial state vector and result
 */
void DEInteg(function<void(double, double **, double **)> func, double t,
             double tout, double relerr, double abserr, int n_eqn, double **y);
