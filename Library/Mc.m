function f=Mc(z,m,q,type,ntrms)
% f=Mc(z,m,q,type,ntrms)

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

if nargin<5, ntrms=50; end
if nargin<4, type=1; end
z=z(:); nz=length(z); q2=sqrt(q);
k=0:ntrms-1; sgn=cos(pi*k)';
mh=1+fix(m/2); u1=q2*exp(-z); u2=q2*exp(z);  
if mod(m,2)==0 % even index   
  r=m/2; [a,c]=matue(q,1,1,ntrms); 
  u=c(:,mh); [ubig,s]=max(abs(u));
  p0=(-1)^r/u(s)*sgn;
  if s==1, p0=p0/2; end; s=s-1; 
  %f=(besselj(k-s,u1).*besselj(k+s,u2)+...
  %   besselj(k+s,u1).*besselj(k-s,u2))*...
  %   (p0.*u);
 f=(besselj(k-s,u1).*Z(k+s,type,u2)+...
     besselj(k+s,u1).*Z(k-s,type,u2))*...
     (p0.*u);
     
else           % odd index
  [a,c]=matue(q,2,1,ntrms); r1=(m+1)/2;
  u=c(:,mh); [ubig,s]=max(abs(u));
  p0=-(-1)^r1/u(s)*sgn; s=s-1;
  %f=(besselj(k-s,u1).*besselj(k+s+1,u2)+...
  %   besselj(k+s+1,u1).*besselj(k-s,u2))*...
  %   (p0.*u); 
  f=(besselj(k-s,u1).*Z(k+s+1,type,u2)+...
     besselj(k+s+1,u1).*Z(k-s,type,u2))*...
     (p0.*u); 
end

%========================================

function f=bes(n,z)
% f=bes(n,z) calls function besselj
% for vector values of n and z.
% f(i,j)= besselj(n(j),z(i))
f= besselj(ones(size(z(:)))*n(:)',z(:)*ones(size(n(:)')));

%========================================

function f=besy(n,z)
% f=besy(n,z) calls function besselj
% for vector values of n and z.
% f(i,j)= bessely(n(j),z(i))
f= bessely(ones(size(z(:)))*n(:)',z(:)*ones(size(n(:)')));

%=========================================

function f=besh(n,k,z)
% f=besh(n,k,z) calls function besselj
% for vector values of n and z.
% f(i,j)= besselj(n(j),k,z(i))
f= besselh(ones(size(z(:)))*n(:)',k,z(:)*ones(size(n(:)')));

%=========================================


function v=Z(n,type,x)
% v=Z(n,type,x)
switch type
  case 1, v=bes(n,x);
  case 2, v=besy(n,x);  
  case 3, v=besh(n,1,x);
  case 4, v=besh(n,2,x);
end