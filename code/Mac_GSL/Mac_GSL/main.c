#include <stdio.h>
#include <math.h>

#include <gsl/gsl_sf_bessel.h>
//#include <mtc_carlo_uni.c>



 int main (void)
{
    double x, y;
    x = 5.0;
    y = gsl_sf_bessel_J0 (x);
    printf ("J0(%g) = %.18e\n", x, y);
    printf ("x:%g,\n y:%.18e\n", x, y);
    return 0;
}

