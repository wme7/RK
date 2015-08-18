function [alpha,beta]=butcher2shuosher(A,b,r);
%function [alpha,beta]=butcher2shuosher(A,b,r);
%
%By David Ketcheson
%
%Generate Shu-Osher form of an explicit Runge-Kutta method,
%given its Butcher form and radius of absolute monotonicity
%
%For an m-stage method, A should be an m x m matrix
%and b should be a column vector of length m.
%
%Note that MATLAB indexes from 1, while the Shu-Osher coefficients
%are usually indexed from zero.

if nargin<3 r = am_radius(A,b); end

s=size(A,1);
K=[A;b'];
G=eye(s)+r*A;
beta=K/G;
alpha=r*beta;
for i=2:s+1
  alpha(i,1)=1-sum(alpha(i,2:s));
end
