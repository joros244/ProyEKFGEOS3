//
// Created by joros on 27/04/2023.
//

#include "../include/mjday.h"
#include <cmath>

double mjday(int yr, int mon, int day, int hr, int min, double sec) {
    double jd =
            367.0 * yr - floor((7 * (yr + floor((mon + 9) / 12.0))) * 0.25) + floor(275 * mon / 9.0) + day + 1721013.5 +
            ((sec / 60.0 + min) / 60.0 + hr) / 24.0;
    return jd - 2400000.5;
}
