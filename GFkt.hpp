#ifndef GFKT_HPP
#define GFKT_HPP
typedef double (*FunctionPointer)(double,double);
class GFkt
{
  double *u;
  Domain *grid;
  int m,n,dimension;
public:
  GFkt(Domain *);
  GFkt(const GFkt&);
  GFkt(GFkt&&) noexcept;
  ~GFkt();

  double& value(int,int);
  GFkt& operator=(const GFkt&);
  GFkt operator+(const GFkt&) const;
  GFkt operator*(const double) const;
  GFkt dx();
  GFkt dy();
  GFkt Laplacian();
  void show();
  GFkt discretize(FunctionPointer f);
  void Output(const char*);
};
#endif
