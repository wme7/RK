function p=rk_order(A,b,c)
%function p=rk_order(A,b,c)
%
%By David Ketcheson
%
%Determine order of a RK method, up to sixth order
%Order conditions from text of Hairer, Norsett, & Wanner
%
%For an m-stage method, input A should be a m x m matrix; 
%b and c should be column vectors of length m

eps=1.e-14;
m=length(b); % # of stages
em=ones(m,1);
p=0;

if sum(b)-1<eps
  p=1;
end

z(1)=sum(A'*b)-1/2;
if (p==1 && abs(z(1))<eps) p=2; end

z(1)=c'.^2*b-1/3;
z(2)=b'*A^2*em-1/6;
if(max(abs(z))<eps && p==2) p=3; end

z(1)=b'*c.^3-1/4;
z(2)=(b.*c)'*A^2*ones(m,1)-1/8;
z(3)=b'*A*c.^2-1/12;
z(4)=b'*A^2*c-1/24;
if(max(abs(z))<eps && p==3) p=4; end

z(1)=c'.^4*b-1/5;
z(2)=(b.*c.^2)'*A*c-1/10;
z(3)=b'*(A*c).^2-1/20;
z(4)=(b.*c)'*A*c.^2-1/15;
z(5)=b'*A*c.^3-1/20;
z(6)=(b.*c)'*A^2*c-1/30;
z(7)=b'*A*diag(c)*A*c-1/40;
z(8)=b'*A^2*c.^2-1/60;
z(9)=b'*A^3*c-1/120;
if(max(abs(z))<eps && p==4) p=5; end

if p==5 
  z(1)=c'.^5*b-1/6;
  z(2)=b'*diag(c).^3*A*c-1/12;
  z(3)=b'*diag(c)*(A*c).^2-1/24;
  z(4)=b'*diag(c).^2*A*c.^2-1/18;
  z(5)=b'*((A*c.^2).*(A*c))-1/36;
  z(6)=b'*diag(c)*A*c.^3-1/24;
  z(7)=b'*A*c.^4-1/30;
  z(8)=b'*diag(c).^2*A^2*c-1/36;
  z(9)=b'*((A^2*c).*(A*c))-1/72;
  z(10)=b'*diag(c)*A*diag(c)*A*c-1/48;
  z(11)=b'*A*diag(c).^2*A*c-1/60;
  z(12)=b'*A*(A*c).^2-1/120;
  z(13)=b'*diag(c)*A^2*c.^2-1/72;
  z(14)=b'*A*diag(c)*A*c.^2-1/90;
  z(15)=b'*A^2*c.^3-1/120;
  z(16)=b'*diag(c)*A^3*c-1/144;
  z(17)=b'*A*diag(c)*A^2*c-1/180;
  z(18)=b'*A^2*diag(c)*A*c-1/240;
  z(19)=b'*A^3*c.^2-1/360;
  z(20)=b'*A^4*c-1/720;
  if(max(abs(z))<eps) p=6; print('This method has order at least six'); end
end
