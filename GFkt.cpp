#include "Curvebase.hpp"
#include "Boundary_curves.hpp"
#include "Domain.hpp"
#include "GFkt.hpp"
#include <iostream>
//p +1 or -1
#define boundary_y(i,j,p)						\
  x_eta=(3*grid->X(i,j)-4*grid->X(i,j+p)+grid->X(i,j+2*p))*(n+1)/3;	\
  y_eta=(3*grid->Y(i,j)-4*grid->Y(i,j+p)+grid->Y(i,j+2*p))*(n+1)/3;	\
  u_eta=(3*  value(i,j)-4*  value(i,j+p)+  value(i,j+2*p))*(n+1)/3;	
#define boundary_x(i,j,p)						\
  x_ksi=(3*grid->X(i,j)-4*grid->X(i+p,j)+grid->X(i+2*p,j))*(m+1)/3;\
  y_ksi=(3*grid->Y(i,j)-4*grid->Y(i+p,j)+grid->Y(i+2*p,j))*(m+1)/3;\
  u_ksi=(3*  value(i,j)-4*  value(i+p,j)+  value(i+2*p,j))*(m+1)/3;
#define inner_eta(i)\
  x_eta=( grid->X(i,j+1)-grid->X(i,j-1) )*(n+1)/2;\
  y_eta=( grid->Y(i,j+1)-grid->Y(i,j-1) )*(n+1)/2;\
  u_eta=(   value(i,j+1)-  value(i,j-1) )*(n+1)/2;
#define inner_ksi(j)\
  x_ksi=( grid->X(i+1,j)-grid->X(i-1,j) )*(m+1)/2;\
  y_ksi=( grid->Y(i+1,j)-grid->Y(i-1,j) )*(m+1)/2;\
  u_ksi=(   value(i+1,j)-  value(i-1,j) )*(m+1)/2;
#define derivative_x(i,j)\
  detJ=x_ksi * y_eta - y_ksi * x_eta;\
  tmp.u[i+j*(m+1)]=1/detJ*(u_ksi*y_eta-u_eta*y_ksi);
#define derivative_y(i,j)\
  detJ=x_ksi * y_eta - y_ksi * x_eta;\
  tmp.u[i+j*(m+1)]=1/detJ*(-u_ksi*x_eta+u_eta*x_ksi);
  
  
  

GFkt::GFkt(Domain *grid_):m(0),n(0),dimension(0),u(nullptr),grid(grid_)
{
  m=grid->xsize();//number of segments
  n=grid->ysize();
  dimension=(m+1)*(n+1);
  u=new double[dimension];
  std::fill(u,u+dimension,0.0);
}

GFkt::GFkt(const GFkt& U):u(nullptr),grid(U.grid)
{
  u=new double[dimension];
  std::copy(U.u,U.u+dimension,u);
}

GFkt::GFkt(GFkt&& B) noexcept:u(B.u),grid(B.grid)
{
  B.grid=nullptr;
  B.u=nullptr;
}

GFkt::~GFkt()
{
  if(u!=nullptr) delete [] u;
}

double& GFkt::value(int ix,int iy)
{
  return u[ix+iy*(m+1)];
}

GFkt& GFkt::operator=(const GFkt& B)
{
  if(this!=&B)
    {
      if(grid!=B.grid)
	{
	  std::cerr<<"different domain for assignment";
	  exit(EXIT_FAILURE);
	}
      std::copy(B.u,B.u+dimension,u);
    }
  return *this;
}

GFkt GFkt::operator+(const GFkt& B) const
{
  if(grid==B.grid)
    {
      GFkt tmp(grid);
      
      for(int i=0;i<(m+1)*(n+1);i++)
	{
	  tmp.u[i]=u[i]+B.u[i];
	}
      return std::move(tmp);
    }
  else
    {
      std::cerr<<"They have different domain!";
      exit(EXIT_FAILURE);
    }
}

GFkt GFkt::operator*(const double a) const
{
  GFkt tmp(grid);
  for(int i=0;i<(m+1)*(n+1);i++)
    {
      tmp.u[i]=u[i]*a;
    }
  return std::move(tmp);
}

GFkt operator*(double a,const GFkt& Q)
{
  GFkt tmp=Q*a;
  return std::move(tmp);
}
GFkt GFkt::dx()
{
  GFkt tmp(grid);//temporary for result
  double x_ksi,y_ksi,u_ksi,x_eta,y_eta,u_eta,detJ;

  for(int j=1;j<n;j++)
    {
      //inner points
      for(int i=1;i<m;i++)
      {
	inner_ksi(j);
	inner_eta(i);
	derivative_x(i,j);
      }
      //i=0
      boundary_x(0,j,+1);
      inner_eta(0);
      derivative_x(0,j);
      //i=m i.e. largest
      boundary_x(m,j,-1);
      inner_eta(m);
      derivative_x(m,j)
    }
  for(int i=1;i<m;i++)
    {
      //j=0
      boundary_y(i,0,+1);
      inner_ksi(0);
      derivative_x(i,0);
      //j=n i.e.largest
      boundary_y(i,n,-1);
      inner_ksi(n);
      derivative_x(i,n);
    }
  //i=0,j=0
  boundary_x(0,0,+1);
  boundary_y(0,0,+1);
  derivative_x(0,0)
  //i=0,j=n
  boundary_x(0,n,+1);
  boundary_y(0,n,-1);
  derivative_x(0,n);
  //i=m,j=0
  boundary_x(m,0,-1);
  boundary_y(m,0,+1);
  derivative_x(m,0);
  //i=m,j=n
  boundary_x(m,n,-1);
  boundary_y(m,n,-1);
  derivative_x(m,n);
  return std::move(tmp);
}

GFkt GFkt::dy()
{
  GFkt tmp(grid);//temporary for result
  double x_ksi,y_ksi,u_ksi,x_eta,y_eta,u_eta,detJ;
  for(int j=1;j<n;j++)
    {
      //inner points
      for(int i=1;i<m;i++)
      {
	inner_ksi(j);
	inner_eta(i);
	derivative_y(i,j);
      }
      //i=0
      boundary_x(0,j,+1);
      inner_eta(0);
      derivative_y(0,j);
      //i=m i.e. largest
      boundary_x(m,j,-1);
      inner_eta(m);
      derivative_y(m,j)
    }
  for(int i=1;i<m;i++)
    {
      //j=0
      boundary_y(i,0,+1);
      inner_ksi(0);
      derivative_y(i,0);
      //j=n i.e.largest
      boundary_y(i,n,-1);
      inner_ksi(n);
      derivative_y(i,n);
    }
  //i=0,j=0
  boundary_x(0,0,+1);
  boundary_y(0,0,+1);
  derivative_y(0,0)
  //i=0,j=n
  boundary_x(0,n,+1);
  boundary_y(0,n,-1);
  derivative_y(0,n);
  //i=m,j=0
  boundary_x(m,0,-1);
  boundary_y(m,0,+1);
  derivative_y(m,0);
  //i=m,j=n
  boundary_x(m,n,-1);
  boundary_y(m,n,-1);
  derivative_y(m,n);
  return std::move(tmp);
}

GFkt GFkt::Laplacian()
{
  GFkt ux(grid),uxx(grid),uy(grid),uyy(grid),tmp(grid);
  ux=dx();
  uy=dy();
  uxx=ux.dx();
  uyy=uy.dy();
  tmp=uxx+uyy;

  return std::move(tmp);
}

void GFkt::show()
{

  for(int j=0;j<n+1;j++)//here it's y index
    {
      for(int i=0;i<m+1;i++)
	std::cout<<u[i+(n-j)*(m+1)]<<"\t";
      std::cout<<"\n";
    }
}

GFkt GFkt::discretize(FunctionPointer f) const
{
  GFkt tmp(grid);
  for(int i=0;i<m+1;i++)
    for(int j=0;j<n+1;j++)
      {
	tmp.u[i+j*(m+1)]=f(grid->X(i,j),grid->Y(i,j));
      }
  return std::move(tmp);
}

void GFkt::Output(const char *name)
{
  FILE *fp;
  fp=fopen(name,"wb");
  fwrite(u,sizeof(double),(m+1)*(n+1),fp);
  fclose(fp);
}
