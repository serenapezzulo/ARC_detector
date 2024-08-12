#ifndef ArcQuarticEquation_H
#define ArcQuarticEquation_H 1

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <math.h>


class ArcQuarticEquation{

 public:
  ArcQuarticEquation();
  virtual ~ArcQuarticEquation();


  int gsl_poly_complex_solve_quartic (double a, double b, double c, double d,
				      gsl_complex* z0, gsl_complex* z1,
				      gsl_complex* z2, gsl_complex* z3);

 private:
  
};
#endif
