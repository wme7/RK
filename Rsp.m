function [gamma,R]=Rsp(s,p)
%function [gamma,R]=Rsp(s,p)
%
%By David Ketcheson
%
%Returns the optimal absolutely monotonic polynomial of degree s 
%and order of accuracy p
%gamma contains the coefficients of the Taylor series about z=-r
%To construct the polynomial, use:
%>  syms z phi
%>  phi=simplify(sum((1.+z/R).^(0:s).*gamma));
%
%Uses the MATLAB optimization toolbox
    
%=========================================================
%Set options for linprog
opts=optimset('TolX',1.e-15,'TolFun',1.e-15,'MaxIter',10000000,...
               'LargeScale','on','Simplex','off','Display','off');
acc=1.e-15; %Accuracy of bisection search
%=========================================================

if p==s %In this case, the optimal polynomial is just the Taylor polynomial
  R=1;
  for i=0:p
    d(i+1)=R^i;
    for j=0:s
      B(i+1,j+1)=prod(j-(0:i-1));
    end
  end
  gamma=(B\ones(s+1,1))';
 
else
  M=s+1;
  rmax=s-p+1.0001;
  rmin=0;
  r=rmax;  %Initial guess
  c=zeros(M,1);
  clear B d;
  
  while (rmax-rmin>acc) %Find R by bisection
    %Set up and improve conditioning of equality constraints
    for i=0:p
      rescale=r^i; d(i+1)=r^i/rescale;
      for j=0:s
        B(i+1,j+1)=prod(j-(0:i-1))/rescale;
      end
    end
    %Test feasibility for this value of r
    [x,lambda,exitflag]=linprog(c,[],[],B,d,zeros(M,1),zeros(M,1)+1.e6,c,opts);
    if exitflag==1;
      rmin=r; r=(r+rmax)/2;
    else
      rmax=r; r=(rmin+r)/2;
    end
  end

  %Now get a feasible solution so we have the coefficients of the method
  R=rmin;
  for i=0:p
    rescale=R^i;
    d(i+1)=R^i/rescale;
    for j=0:s
      B(i+1,j+1)=prod(j-(0:i-1))/rescale;
    end
  end
  [gamma,lambda,exitflag]=linprog(c,[],[],B,d,zeros(M,1),zeros(M,1)+1.e6,c,opts);
end
