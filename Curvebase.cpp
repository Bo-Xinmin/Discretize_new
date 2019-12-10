#include <iostream>
#include <cstdlib>
#include <cmath>
#include "Curvebase.hpp"


double Curvebase::unitlength(double p)
{
  return std::sqrt(dxp(p)*dxp(p)+dyp(p)*dyp(p));
}
double Curvebase::I1(double a,double b)
{
  double g = 0.5*(a+b);
  double I = (b-a)*(unitlength(a)+4*unitlength(g)+unitlength(b))/6;
  return I;
}
double Curvebase::asi(double a, double b,double tol)
{
  double g = 0.5*(a+b);
  double I2 = I1(a,g)+I1(g,b);
  double errest = fabs(I1(a,b)-I2);
  if (errest < 15*tol)
    {
      return I2;
    }
  return asi(a,g,0.5*tol) + asi(g,b,0.5*tol);
}
double Curvebase::integrate(double p)
{
  return asi(pmin,p);
}

// CONSTRUCTOR
Curvebase::Curvebase(double a,double b,bool r):
  pmin(a), pmax(b), rev(r){}

double Curvebase::newton(double s)
{
  double p_i = 0.5*(pmax+pmin); // Initial guess
  // Initial values of function and derivative
  double f = integrate(p_i) - s*length;
  double df = unitlength(p_i);
  while (fabs(f) > TOL)
    {
      p_i -= f/df;// Update p_i
      // Calulate new function and derivative with new p_i
      f = integrate(p_i) - s*length;
      df = unitlength(p_i);
    }
  return p_i;
}

double Curvebase::x(double s)
{
  //if(length==0)
  //  length = integrate(pmax);
  if(s<0 || s>1) // Check if s in not in the interval
    {
      std::cerr<<"the distribution parameter can only be between 0 and 1";
      exit(EXIT_FAILURE);
    }
  double p = newton(s);
  return xp(p);
}
double Curvebase::y(double s)
{
  if(s<0 || s>1) // Check if s in not in the interval
    {
      std::cerr<<"the distribution parameter can only be between 0 and 1";
      exit(EXIT_FAILURE);
    }
  double p = newton(s);
  return yp(p);
}

