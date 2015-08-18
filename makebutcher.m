function [A,b,c,r]=makebutcher(name,s)
%function [A,b,c,r]=makebutcher(name,s)
%
%By David Ketcheson
%
%Set up Butcher arrays A,b,c for various methods
%Also returns SSP coefficient r
%For families of methods, optional input s is the number of stages

if nargin<2 s=1; end

switch name
%=================SSP Methods=========================

  %=================Explicit Methods=========================
  case 'FE11'
    %Forward Euler
    s=1; r=1;
    A=[0];
    b=[1]'; c=[0]';
  
  case 'SSP22'
    s=2; r=1;
    A=[0 0; 1 0];
    b=[1/2 1/2]'; c=sum(A,2);

  case 'SSP42'
    s=4; r=3;
    A=[0 0 0 0; 1/3 0 0 0; 1/3 1/3 0 0; 1/3 1/3 1/3 0];
    b=1/4*ones(m,1); c=sum(A,2);

  case 'SSP33'
    s=3; r=1;
    A=[0 0 0; 1 0 0; 1/4 1/4 0];
    b=[1/6 1/6 2/3]'; c=sum(A,2);

  case 'SSP43'
    s=4; r=2;
    A=[0 0 0 0; 1/2 0 0 0; 1/2 1/2 0 0; 1/6 1/6 1/6 0];
    b=[1/6 1/6 1/6 1/2]'; c=sum(A,2);

  case 'SSP104'
    s=10; r=6;
    alpha0=diag(ones(1,s-1),-1);
    alpha0(6,5)=2/5; alpha0(6,1)=3/5;
    beta0 =1/6*diag(ones(1,s-1),-1);
    beta0(6,5)=1/15;
    A=(eye(s)-alpha0)\beta0;
    b=1/10*ones(s,1); c=sum(A,2);

  case 'rSSPs2'
    %Rational (optimal, low-storage) s-stage 2nd order SSP
    if s<2 error('Explicit second order SSP family requires s>=2'); end
    r=s-1;
    alpha=[zeros(1,s);eye(s);];
    alpha(s+1,s)=(s-1)/s;
    beta=alpha/r;
    alpha(s+1,1)=1/s;
    A=(eye(s)-alpha(1:s,:))\beta(1:s,:);
    b=beta(s+1,:)+alpha(s+1,:)*A; b=b';
    c=sum(A,2);

  case 'rSSPs3'
    %Rational (optimal, low-storage) s^2-stage 3rd order SSP
    if round(sqrt(s))~=sqrt(s) || s<4
       error('Explicit third order SSP family requires s=n^2, n>1'); 
    end
    n=s^2; r=n-s;
    alpha=[zeros(1,n);eye(n);];
    alpha(s*(s+1)/2+1,s*(s+1)/2)=(s-1)/(2*s-1);
    beta=alpha/r;
    alpha(s*(s+1)/2+1,(s-1)*(s-2)/2+1)=s/(2*s-1);
    A=(eye(n)-alpha(1:n,:))\beta(1:n,:);
    b=beta(n+1,:)+alpha(n+1,:)*A; b=b';
    c=sum(A,2);


  %=================Implicit Methods=========================
  case 'BE11'
    %Backward Euler
    s=1; r=1.e10;
    A=[1];
    b=[1]'; c=[1]';
  
  case 'SDIRK34' %3-stage, 4th order singly diagonally implicit (SSP)
    s=3; r=1.7588;
    g=0.5*(1-cos(pi/18)/sqrt(3)-sin(pi/18));
    q=(0.5-g)^2;
    A=[g     0    0
       0.5-g g    0
       2*g  1-4*g g];
    b=[1/(24*q) 1-1/(12*q) 1/(24*q)]';
    c=sum(A,2);

   case 'ISSPm2'
    %Optimal DIRK SSP schemes of order 2
    r=2*s;
    i=repmat((1:s)',1,s); j=repmat(1:s,s,1);
    A=1/s*(j<i) + 1/(2*s)*(i==j);
    b=1/s*ones(s,1);
    c=sum(A,2);

  case 'ISSPs3'
    %Optimal DIRK SSP schemes of order 3
    if s<2 error('Implicit third order SSP schemes require s>=2'); end
    r=s-1+sqrt(s^2-1);
    i=repmat((1:s)',1,s); j=repmat(1:s,s,1);
    A=1/sqrt(s^2-1)*(j<i) + 0.5*(1-sqrt((s-1)/(s+1)))*(i==j);
    b=1/s*ones(s,1);
    c=sum(A,2);


%=================Classical Methods=========================

  %Gauss-Legendre methods -- order 2s
  case 'GL1'
    r=2; A=1/2; b=1; c=1/2;
  case 'GL2'
    r=0;
    A=[1/4 1/4-sqrt(3)/6
       1/4+sqrt(3)/6 1/4];
    b=[1/2 1/2]';
    c=[1/2-sqrt(3)/6 1/2+sqrt(3)/6]';
  case 'GL3'
    r=0;
    A=[5/36 (80-24*sqrt(15))/360 (50-12*sqrt(15))/360
       (50+15*sqrt(15))/360 2/9  (50-15*sqrt(15))/360
       (50+12*sqrt(15))/360 (80+24*sqrt(15))/360 5/36];
    b=[5/18 4/9 5/18]';
    c=[(5-sqrt(15))/10 1/2 (5+sqrt(15))/10]';

  %Radau IA methods -- order 2s-1
  case 'RIA1'
    r=1;
    A=1; b=1; c=0;
  case 'RIA2'
    r=0;
    A=[1/4 -1/4
       1/4 5/12];
    b=[1/4 3/4]';
    c=[0 2/3]';
  case 'RIA3'
    r=0;
    A=[1/9 (-1-sqrt(6))/18 (-1+sqrt(6))/18
       1/9 (88+7*sqrt(6))/360 (88-43*sqrt(6))/360
       1/9 (88+43*sqrt(6))/360 (88-7*sqrt(6))/360];
    b=[1/9 (16+sqrt(6))/36 (16-sqrt(6))/36]';
    c=[0 (6-sqrt(6))/10 (6+sqrt(6))/10]';

  %Radau IIA methods -- order 2s-1
  case 'RIIA1'
    r=1;
    A=1; b=1; c=1;
  case 'RIIA2'
    r=0;
    A=[5/12 -1/12
       3/4 1/4];
    b=[3/4 1/4]';
    c=[1/3 1]';
  case 'RIIA3'
    r=0;
    A=[(88-7*sqrt(6))/360 (296-169*sqrt(6))/1800 (-2+3*sqrt(6))/225
       (296+169*sqrt(6))/1800 (88+7*sqrt(6))/360 (-2-3*sqrt(6))/225
       (16-sqrt(6))/36 (16+sqrt(6))/36 1/9];
    b=[(16-sqrt(6))/36 (16+sqrt(6))/36 1/9 ]';
    c=[(4-sqrt(6))/10 (4+sqrt(6))/10 1];

  %Lobatto IIIA methods -- order 2s-2
  case 'LIIIA2'
    r=0;
    A=[0 0
       1/2 1/2];
    b=[1/2 1/2]';
    c=[0 1]';
  case 'LIIIA3'
    r=0;
    A=[0 0 0
       5/24 1/3 -1/24
       1/6 2/3 1/6];
    b=[1/6 2/3 1/6]';
    c=[0 12 1];

%===================Miscellaneous Methods================

  case 'Mid22'
    %Midpoint 22 method
    s=2; r=0.5;
    A=[0  0
       1/2 0];
    b=[0 1]'; c=[0 1/2]';

  case 'MTE22'
    %Minimal truncation error 22 method (Heun)
    s=2; r=0.5;
    A=[0  0
       2/3 0];
    b=[1/4 3/4]'; c=[0 2/3]';

  case 'CN22'
    %Crank-Nicholson
    s=2; r=2;
    A=[0  0
       1/2 1/2];
    b=[1/2 1/2]'; c=[0 1]';

  case 'Heun33'
    s=3; r=0;
    A=[0 0 0; 1/3 0 0; 0 2/3 0];
    b=[1/4 0 3/4]'; c=sum(A,2);

  case 'RK44'  %Classical fourth order
    s=4; r=0;
    A=[0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
    b=[1/6 1/3 1/3 1/6]'; c=sum(A,2);

%====================DSRK Methods========================

  case 'DSso2'
    %CBM's DSRKso2
    s=2; isdsrk=1;
    A=[3/4 -1/4
        1    0];
    W=[1/2   0
       1     0];
    b=[1 0]'; c=[1/2 1]';

  case 'DSRK2'
    %CBM's DSRK2
    s=2; isdsrk=1;
    A=[1/2 -1/2
       1/2  1/2];
    W=[ 0    0
       1/2  1/2];
    b=[1/2 1/2]'; c=[0 1]';

  case 'DSRK3'
    %Zennaro's DSRK3
    s=3; isdsrk=1;
    A=[5/2 -2 -1/2
       -1  2  -1/2
       1/6 2/3 1/6];
    W=[ 0  0  0
       7/24 1/6 1/24
       1/6 2/3 1/6];
    b=[1/6 2/3 1/6]'; c=[0 1/2 1]';
 %===================="Non-SSP" Methods of Wong & Spiteri========================
  case 'NSSP21'
    m=2; r=0;
    A=[0   0
       3/4 0];
    b=[0 1]'; c=[0 3/4]';

  case 'NSSP32'
    m=3; r=0;
    A=[0   0 0
       1/3 0 0
       0   1 0];
    b=[1/2 0 1/2]'; c=[0 1/3 1]';

  case 'NSSP33'
    m=3; r=0;
    A=[0    0    0
       -4/9 0    0
       7/6  -1/2 0];
    b=[1/4 0 3/4]'; c=[0 -4/9 2/3]';

  case 'NSSP53'
    m=5; r=0;
    A=[0 0 0 0 0
       1/7 0 0 0 0
       0 3/16 0 0 0
       0 0 1/3 0 0
       0 0 0 2/3 0];
    b=[1/4 0 0 0 3/4]'; c=[0 1/7 3/16 1/3 2/3]';

   
end
