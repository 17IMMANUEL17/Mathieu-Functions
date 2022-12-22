function f=Mc3big(z,m,q,ntrms)
% f=Mc3big(z,m,q,ntrms)

% Modified Mathieu function Mc(z,m,q,type)
% defined by 20.4.7, 20.6.7 & 20.6.8 of
% 'Handbook of Mathematical Functions'
% by Milton Abramowitz and Irene Stegun
% z     - vector of argument values
% m     - integer function order
% q     - scalar parameter value
% type  - 1,2,3, or 4 to use besselj,
%         bessely,hankel1, or hankel2
% ntrms - number of series terms used
% f     - vector of function values

%          by H.W. 3/20/06

if nargin<4, ntrms=50; end
z=z(:); nz=length(z); q2=sqrt(q);
k=(0:ntrms-1); sgn=cos(pi*k)';
mh=1+fix(m/2); uu=sqrt(4*q)*cosh(z); %u1=q2*exp(-z); u2=q2*exp(z);  
if mod(m,2)==0 % even index   
  r=m/2; [a,c]=matue(q,1,1,ntrms); 
  kk=2*(0:ntrms-1);
  u=c(:,mh); %[ubig,s]=max(abs(u));
  %p0=(-1)^r/u(s)*sgn;
  %p0=(-1)^r*sgn;
  %if s==1, p0=p0/2; end; s=s-1; 
  %f=(besselj(k-s,u1).*besselj(k+s,u2)+...
  %   besselj(k+s,u1).*besselj(k-s,u2))*...
  %   (p0.*u);
 %f=(besselj(k-s,u1).*Z(k+s,type,u2)+...
 %    besselj(k+s,u1).*Z(k-s,type,u2))*...
 %    (p0.*u);
 f=besselh(kk,1,uu)*(sgn.*u);
     
else           % odd index
  [a,c]=matue(q,2,1,ntrms); r1=(m+1)/2;
  u=c(:,mh); %[ubig,s]=max(abs(u));
  %p0=-(-1)^r1/u(s)*sgn; s=s-1;
  %f=(besselj(k-s,u1).*besselj(k+s+1,u2)+...
  %   besselj(k+s+1,u1).*besselj(k-s,u2))*...
  %   (p0.*u); 
  %f=(besselj(k-s,u1).*Z(k+s+1,type,u2)+...
  %   besselj(k+s+1,u1).*Z(k-s,type,u2))*...
  %   (p0.*u); 
  kk=1+2*(0:ntrms-1);
  f=besselh(kk,1,uu)*(sgn.*u);
end

%======================================

% function v=Z(n,type,x)
% % v=Z(n,type,x)
% switch type
%   case 1, v=bes(n,x);
%   case 2, v=besy(n,x);  
%   case 3, v=besh(n,1,x);
%   case 4, v=besh(n,2,x);
% end

%======================================

function [a,c,v]=matue(q,order,ftype,ntrms)
% [a,c,v]=matue(q,order,ftype,ntrms)

% For a given value of q, this function computes
% values of parameter a in the Mathieu equation 
%     y"(x)+(a-2*q*cos(2*x))*y(x)=0
% The coefficients in the Fourier series
% expansions of the Mathieu functions
% ce(x,n,q) and se(x,n,q) are also obtained.
% q     - parameter in the Mathieu equation
%         proportional to the square of the
%         related natural frequency of an
%         elliptic membrane. This value can
%         be complex.
% order - 1 for even index, 2 for odd index
% ftype - 1 for even valued, 2 for odd-valued
% ntrms - number of terms used to approximate
%         the Fourier expansion of ce or se
% a     - eigenvalues which produce periodic
%         solutions of the Mathieu equation
% c     - matrix of eigenvectors defining the
%         Fourier coefficients of the cosine
%         or sine series
% v     - tridiagonal matrix which determines 
%         the eigenvalues and eigenvectors

%         HBW 12/01/04

if nargin<4, ntrms=50; end
if nargin==0, q=10; order=1; ftype=1; end
v=diag(q*ones(ntrms-1,1),1); v=v+v.';
if ftype==1;  % even valued
  if order==1 % even order
    v=v+diag((2*(0:ntrms-1)).^2);
    v(1,2)=sqrt(2)*q; v(2,1)=v(1,2);     
  else        % odd order
    v=v+diag((1+2*(0:ntrms-1)).^2);
    v(1,1)=1+q;
  end
else % ftype==2, odd valued
  if order==1  % even order
    v=v+diag((2*(1:ntrms)).^2);
  else         % odd order 
    v=v+diag((1+2*(0:ntrms-1)).^2);  
    v(1,1)=1-q;
  end
end
if abs(q)>1000,v=v/q; end
[c,a]=eig(v); 
if abs(q)>1000, a=a*q; end
a=diag(a); [dumy,k]=sort(real(a));
a=a(k); c=c(:,k);
if order==1 & ftype==1
  c(1,:)=c(1,:)/sqrt(2);
end    
% Normalize the series coefficients
cv=sum(conj(c).*c);
if order==1 & ftype==1
  cv=cv+conj(c(1,:)).*c(1,:);
end
if~isreal(q) % Complex eigenvector case
  c=ones(ntrms,1)*(1./sqrt(cv)).*c;
  return
end
% Case where q is real. Make ce(0,q)
% or se'(0,q) positive.
if ftype==1 % even valued
  sgn=sign(sum(c)); sgn=sgn+(sgn==0);  
  cv=sgn./sqrt(cv);
else        % odd valued
  if order==1, n=2*(1:ntrms);
  else, n=2*(1:ntrms)-1; end
  sgn=sign(n*c); sgn=sgn+(sgn==0); 
  cv=(sgn./sqrt(cv));
end  
c=c.*cv(ones(ntrms,1),:); 
