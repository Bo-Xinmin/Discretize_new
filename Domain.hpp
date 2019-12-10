#ifndef DOMAIN_HPP
#define DOMAIN_HPP
class Domain
{
  Curvebase *sides[4];
  double *x,*y; //Array of coordinates of each inner node
  int m,n; //Number of internal nodes
public:
  Domain(Curvebase&,Curvebase&,Curvebase&,Curvebase&); //Constructor
  Domain(const Domain&); //Copy constructor
  Domain& operator=(const Domain&); //Copy assignment
  ~Domain(); // Destructor
  int xsize(){return m;}
  int ysize(){return n;}
  bool grid_valid(){return (m*n)!=0;}
  void generate_grid(int,int);//number of intervals in each direction
  void Output(); // Writting internal nodes to file
  double X(int,int);//x index,y index
  double Y(int,int);
};
#endif
