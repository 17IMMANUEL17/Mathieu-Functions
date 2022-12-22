function f=Ms(z,m,q,type,ntrms)
% f=Ms(z,m,q,type,ntrms)

% Modified Mathieu function Ms(z,m,q)
% defined by 20.4.7, 20.6.9 & 20.6.10 of
% 'Handbook of Mathematical Functions'
% by Milton Abramowitz and Irene Stegun
% z     - vector of argument values
% m     - integer function order
% q     - scalar parameter value
% ntrms - number of series terms used
% f     - vector of function values

%       by H.W. 4/20/06

if nargin<5, ntrms=50; end
if nargin<4, type=1; end
z=z(:); z1=sqrt(q)*exp(-z); 
z2=sqrt(q)*exp(z); 
if mod(m,2)==0 % even order
  k=1:ntrms; [a,c]=matue(q,1,2,ntrms);
  u=c(:,m/2); r=m/2;  
  [ubig,s]=max(abs(u)); p0=(-1)^r/u(s);
  %f=(bes(k-s,z1).*bes(k+s,z2)-...
  %   bes(k+s,z1).*bes(k-s,z2))*...
  %   (p0*cos(k*pi)'.*u);
  f=(bes(k-s,z1).*Z(k+s,type,z2)-...
     bes(k+s,z1).*Z(k-s,type,z2))*...
     (p0*cos(k*pi)'.*u);
else           % odd order
  k=0:ntrms-1; [a,c,vsave]=matue(q,2,2,ntrms);
  u=c(:,1+fix(m/2)); [ubig,s]=max(abs(u));
  r=(m-1)/2; p0=(-1)^r/u(s); s=s-1;
  %f=(bes(k-s,z1).*bes(k+s+1,z2)-...
  %   bes(k+s+1,z1).*bes(k-s,z2))*...
  %   (p0*cos(k*pi)'.*u);
  f=(bes(k-s,z1).*Z(k+s+1,type,z2)-...
     bes(k+s+1,z1).*Z(k-s,type,z2))*...
     (p0*cos(k*pi)'.*u);
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
