#include "Curvebase.hpp"
#include <cmath>
#ifndef BOUNDARY_CURVES_HPP
#define BOUNDARY_CURVES_HPP
class LowerCurve:public Curvebase
{
  double xp(double p)
  {
    return p;
  }
  double yp(double p)
  {
    if(p<-3)
      return 0.5/(1+exp(-3*(p+6)));
    else
      return 0.5/(1+exp(3*p));
  }
  double dxp(double p)
  {
    return 1.0;
  }
  double dyp(double p)
  {
    if(p<-3)
      {
	double t=-3*(p+6);
	return 1.5*exp(t)/((1+exp(t))*(1+exp(t)));
      }
    else 
      {
	double t=3*p;
	return -1.5*exp(t)/((1+exp(t))*(1+exp(t)));
      }
  }
public:
  LowerCurve(double a,double b):Curvebase(a,b)
  {
    length = integrate(b);
  }
  ~LowerCurve(){};
};

class LeftCurve:public Curvebase
{
  double xp(double p)
  {
    return -10.0;
  }
  double yp(double p)
  {
      return p;
  }
  double dxp(double p)
  {
    return 0.0;
  }
  double dyp(double p)
  {
      return 1.0;
  }
public:
  LeftCurve(double a,double b):Curvebase(a,b)
  {
    length = integrate(b);
  }
  ~LeftCurve(){};
};

class RightCurve:public Curvebase
{
  double xp(double p)
  {
    return 5.0;
  }
  double yp(double p)
  {
      return p;
  }
  double dxp(double p)
  {
    return 0.0;
  }
  double dyp(double p)
  {
      return 1.0;
  }
public:
  RightCurve(double a,double b):Curvebase(a,b)
  {
    length = integrate(b);
  }
  ~RightCurve(){};
};

class UpperCurve:public Curvebase
{
  double xp(double p)
  {
      return p;
  }
  double yp(double p)
  {
    return 3.0;
  }
  double dxp(double p)
  {
      return 1.0;;
  }
  double dyp(double p)
  {
    return 0.0;
  }
public:
  UpperCurve(double a,double b):Curvebase(a,b)
  {
    length = integrate(b);
  }
  ~UpperCurve(){};
};
#endif
