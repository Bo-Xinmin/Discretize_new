#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include "Curvebase.hpp"
#include "Boundary_curves.hpp"
#include "Domain.hpp"
// Constructor
Domain::Domain(Curvebase& s1,Curvebase& s2,Curvebase& s3,Curvebase& s4)
{
  sides[0]=&s1;//Left
  sides[1]=&s2;//Right
  sides[2]=&s3;//Lower
  sides[3]=&s4;//Upper
  m=n=0;
  x=y=nullptr;
}
//Copy constructor
Domain::Domain(const Domain& d):m(d.m),n(d.n),x(nullptr),y(nullptr)
{
  if(m>0 && n>0)
    {
      x=new double[(m+1)*(n+1)];
      y=new double[(m+1)*(n+1)];
      for(int i=0;i<(m+1)*(n+1);i++)
	{
	  x[i]=d.x[i];
	  y[i]=d.y[i];
	}
    }
}
// Copy assignment
Domain& Domain::operator=(const Domain& d)
{
  if(this!=&d)
    {//dont copy to itself
      if(m==d.m && n==d.n)
	for(int i=0;i<(m+1)*(n+1);i++)
	  {
	    x[i]=d.x[i];
	    y[i]=d.y[i];
	  }
      else
	{
	  if(m>0)
	    {
	      delete [] x;
	      delete [] y;
	      x=y=nullptr;
	    }
	  if(m>0)
	    {
	      x=new double[(m+1)*(n+1)];
	      y=new double[(m+1)*(n+1)];
	      for(int i=0;i<(m+1)*(n+1);i++)
		{
		  x[i]=d.x[i];
		  y[i]=d.y[i];
		}
	    }
	}
    }
}
// Destructor
Domain::~Domain()
{
  if(m>0)
    {
      delete [] x;
      delete [] y;
    }
}
void Domain::generate_grid(int m_,int n_)
{
  //m on x;n on y
  //m,n numbers of intervals
  if(m_<=0 || n_<=0)
    {
      std::cerr<<"meaningless number!";
      exit(EXIT_FAILURE);
    }
  else
    {
      if(m>0)
	{
	  delete [] x;
	  delete [] y;
	}
      m=m_;n=n_;
      x=new double[(m+1)*(n+1)];
      y=new double[(m+1)*(n+1)];
      double ksi,eta,sigma;
      double Corner_x[4]={-10.0,5.0,-10.0,5.0};
      //Specialized here, only works in this problem
      double Corner_y[4]={0.0,0.0,3.0,3.0};
      int k;
      const double delta=3.0;
      for(int i=0;i<n+1;i++)//y direction
	{
	for(int j=0;j<m+1;j++)//x direction
	  {
	    ksi=(double)(j)/(double)(m);//ksi,j,m in the x direction
	    
	    //non-uniform distribution
	    sigma=(double)(i)/(double)(n); //eta,i,n in the y direction
	    //eta=1.0+((std::tanh(delta*(sigma-1))))/(std::tanh(delta));
	    eta=sigma;
	    k=i*(m+1)+j;//Lexicograghical rule
	    
	    x[k]=(1-ksi) * (*sides[0]).x(eta) + (ksi) * (*sides[1]).x(eta)
	      + (1-eta) * (*sides[2]).x(ksi) + (eta) * (*sides[3]).x(ksi)
	      - (1-ksi)*(1-eta) * Corner_x[0] - (ksi)*(1-eta) * Corner_x[1]
	      - (1-ksi)*(eta) * Corner_x[2] - (ksi)*(eta) * Corner_x[3];
	    y[k]=(1-ksi) * (*sides[0]).y(eta) + (ksi) * (*sides[1]).y(eta)
	      + (1-eta) * (*sides[2]).y(ksi) + (eta) * (*sides[3]).y(ksi)
	      - (1-ksi)*(1-eta) * Corner_y[0] - (ksi)*(1-eta) * Corner_y[1]
	      - (1-ksi)*(eta) * Corner_y[2] - (ksi)*(eta) * Corner_y[3];
	    // std::cout<<k<<"\t"<<ksi<<"\t"<<eta<<"\t"<<x[k]<<"\t"<<y[k]<<"\n";
	    
	  }
	}
    }
}
void Domain::Output()
{
  FILE *fp;
  //fp=fopen("outfile.bin","wb");
  fp=fopen("Domain.bin","wb");
  fwrite(x,sizeof(double),(m+1)*(n+1),fp);
  fwrite(y,sizeof(double),(m+1)*(n+1),fp);
  fclose(fp);
}
double Domain::X(int ix,int iy)
{
  return x[ix+iy*(m+1)];
}
double Domain::Y(int ix,int iy)
{
  return y[ix+iy*(m+1)];
}
