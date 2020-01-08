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
  ~LowerCurve(){}
};

class StraightLine:public Curvebase
{
  double ax,bx,ay,by;
  double xp(double p)
  {
    return ax*p+bx;
  }
  double yp(double p)
  {
    return ay*p+by;
  }
  double dxp(double p)
  {
    return ax;
  }
  double dyp(double p)
  {
    return ay;
  }
public:
  StraightLine(double ax_,double bx_,double ay_,double by_,double a,double b):Curvebase(a,b),ax(ax_),bx(bx_),ay(ay_),by(by_)
  {
    length = integrate(b);
  }
  ~StraightLine(){}
};
#endif
