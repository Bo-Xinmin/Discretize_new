#ifndef CURVEBASE_HPP
#define CURVEBASE_HPP
//
const double TOL=1e-5;
//
class Curvebase
{
protected:
  double pmin;
  double pmax;
  bool rev;
  double length;
  virtual double xp(double)=0;
  virtual double yp(double)=0;
  virtual double dxp(double)=0;
  virtual double dyp(double)=0;
  double integrate(double);
  double I1(double,double);
  double asi(double, double, double=TOL);
  double unitlength(double);
  double newton(double);
public:
  Curvebase(double, double, bool = true);
  double x(double);
  double y(double);
  virtual ~Curvebase(){};
};

#endif

