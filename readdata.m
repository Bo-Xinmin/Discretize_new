clear;clc;
FILE=fopen("Domain.bin");
d=fread(FILE,'double');
FILE=fopen("dx.bin");
fx=fread(FILE,'double');
FILE=fopen("dy.bin");
fy=fread(FILE,'double');

d=reshape(d,[45,2]);
x=d(:,1);y=d(:,2);
fx_real=0.02*x.*cos(0.01*x.*x).*cos(0.1*x)-0.1*sin(0.1*x).*sin(0.01*x.*x);
diff_x=reshape(fx-fx_real,[9,5]);

