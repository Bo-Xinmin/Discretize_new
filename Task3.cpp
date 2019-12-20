#include "Curvebase.hpp"
#include "Boundary_curves.hpp"
#include "Domain.hpp"
#include "GFkt.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
double func(double x,double y)
{
  return sin(x*0.1*x*0.1)*cos(x*0.1)+y;
}
int main()
{
  LeftCurve L1(0.0,3.0);
  RightCurve L2(0.0,3.0);
  LowerCurve L3(-10.0,5.0);
  UpperCurve L4(-10.0,5.0);

  Domain D(L1,L2,L3,L4);
  D.generate_grid(8,4);//number of inner nodes
  D.Output();

  GFkt u(&D),up(&D),ux(&D),uy(&D),uL(&D);
  u.show();
  
  up=u.discretize(&func);
  up.show();

  ux=up.dx();
  ux.show();
  ux.Output("dx.bin");

  uy=up.dy();
  uy.show();
  uy.Output("dy.bin");

  uL=up.Laplacian();
  uL.show();
  uL.Output("Laplace.bin");
  return 1;
}
