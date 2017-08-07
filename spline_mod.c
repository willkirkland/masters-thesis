#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_spline.h>
/*#include <stdio.h> include for debugging only*/

gsl_spline *spline;
gsl_interp_accel *splineAccel;

void spline_init_(double *xVals, double *yVals, size_t *size)
{
  int status;

  spline = gsl_spline_alloc(gsl_interp_cspline, *size);
  splineAccel = gsl_interp_accel_alloc();
  status = gsl_spline_init(spline, xVals, yVals, *size);
}

void spline_read_(double *inputVal, double *outputVal)
{
  *outputVal = gsl_spline_eval(spline, *inputVal, splineAccel);
}

void spline_free_()
{
  gsl_spline_free(spline);
  gsl_interp_accel_free(splineAccel);
}
