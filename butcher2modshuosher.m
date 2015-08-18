function [lambda,mu1]=butcher2modshuosher(A,b,r);
%function [lambda,mu1]=butcher2modshuosher(A,b,r);
%
%By David Ketcheson
%
%Generate modified Shu-Osher form of a Runge-Kutta method,
%given its Butcher form and radius of absolute monotonicity
%
%For an m-stage method, A should be an m x m matrix
%and b should be a column vector of length m.
%
%Note that MATLAB indexes from 1, while the Shu-Osher coefficients
%are usually indexed from zero.

Aup=triu(A);
if max(abs(Aup))>0 mclass='implicit'; else mclass='explicit'; end

if nargin<3 r = am_radius(A,b); end

s=size(A,1);
K=[A;b'];
G=eye(s)+r*A;
mu1=K/G;
lambda=r*mu1;
for i=1:s+1
  if strcmp(mclass,'implicit') %0 stage is u_n
   alpha(i)=1-sum(lambda(i,1:s));
  end
end

%Eliminate diagonal terms of lambda array (for implicit methods)
lambda=lambda; mu1=mu1;
for i=1:s
  if lambda(i,i)~=1
    fac=1/(1-lambda(i,i));
    mu1(i,:)=fac*mu1(i,:);
    lambda(i,:)=fac*lambda(i,:);
    lambda(i,i)=0.;
  end
end
