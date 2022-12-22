function qplotmatu(func,x,n,q,p)
% A typical call is qplotmatu(@Mc1,atanh(.5), 20,[0,500],1)
% An exponentially plot(for p>=1)is made for any Mathieu
% using a vector of q values and fixed values of x and n.
% Parameter p is used to scale the ordinate values when
% func becomes very small. This is achieved by plotting
% sign(func)*abs(func)^(1/p)
if nargin<5, p=1; end
if nargin==0; func=@Mc1; x=atanh(.5); n=20; q=[0,500]; end
if length(q)==2, q=linspace(q(1),q(2),200); end
f=sq(vecval(func,x,n,q),p);
plot(q,f), grid on
xlabel('q axis'); ylabel('function value')
str=['Plot of ',func2str(func),' for x = ',num2str(x)];
str=[str,' n = ',num2str(n),' and p = ',num2str(p)];
title(str), shg

%================================================

function v=sq(z,p)
% v=sq(z,p)
if nargin<2, p=5; end
z=squeeze(z);
v=sign(z).*abs(z).^(1/p); 

%================================================

function [f,args]=vecval(fname,z,m,q,ntrms)
% [f,args]=vecval(fname,z,m,q,ntrms)
% This function evaluates various Mathieu
% functions allowing vector values for the
% input parameters z,m, and q.
% fname - character string name of a 
%         Mathieu function such as 'Mc2'
% z     - vector of z values
% m     - vector of function orders
% q     - vector of q values
% ntrms - number of series terms
% f       an array in which f(i,j,k) is
%         the derivative of Ms2(z,m,q)
%         evaluated at z(i),m(j),q(k)
% args    a structure returning the
%         input data in elements
if nargin<5, ntrms=50; end
z=z(:); nz=length(z); m=m(:)';
M=length(m); nq=length(q); 
f=zeros(nz,M,nq);
for j=1:M, for k=1:nq
  f(:,j,k)=feval(fname,z,m(j),q(k),ntrms);  
end, end
if nargout==2
  args=struct('z',{z},'m',{m},'q',{q});
end