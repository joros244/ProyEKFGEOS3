#include <stdio.h>
#include <math.h>
#include "include/mjday.h"

using namespace std;

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)


int mjday_01() {

    _assert(fabs(mjday(2023, 4, 27, 18, 17, 9.34) - 60061.76191365718841552734375) < pow(10, -12));

    return 0;
}

int mjday_02() {

    _assert(fabs(mjday(2023, 4, 27) - 60061.0) < pow(10, -12));

    return 0;
}

int mjday_03() {

    _assert(fabs(mjday(0, 0, 0) + 678987.0) < pow(10, -12));

    return 0;
}


int all_tests() {
    _verify(mjday_01);
    _verify(mjday_02);
    _verify(mjday_03);
    return 0;
}


int main() {
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
