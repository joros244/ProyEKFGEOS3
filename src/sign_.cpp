#include "../include/sign_.h"
#include <cmath>

double sign_(double a, double b) {

  if (b >= 0.0) {
    return fabs(a);
  } else {
    return -fabs(a);
  }
}
