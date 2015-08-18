function [r] = cc_radius(A,b)
%function [r] = cc_radius(A,b)
%
%By David Ketcheson
%
%Evaluates the Radius of circle contractivity
%of a Runge-Kutta method, given the Butcher array
%
%For an m-stage method, A should be an m x m matrix
%and b should be a column vector of length m.
%
%Accuracy can be changed by modifying the value of eps.
%Methods with very large radii of a.m. (>1000) will require
%rmax to be increased.

rmax=1000; eps=1.e-13;

if min(b)<=0
  r=0;
else
  m=length(b);
  B=diag(b);
  M=B*A+A'*B-b*b';
  rlo=0; rhi=rmax;
  while rhi-rlo>eps
    r=0.5*(rhi+rlo);
    X=M+B/r;
    if min(eig(X))<-3.e-16
      rhi=r;
    else
      rlo=r;
    end
  end
end
if rhi==rmax % r>=rmax
  error('Error: increase value of rmax in cc_radius.m');
else
  r=rlo;
end
