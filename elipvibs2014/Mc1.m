function f=Mc1(z,m,q,ntrms)
% f=Mc1(z,m,q,ntrms)

% Modified Mathieu function Mc1(z,m,q)
% defined by 20.6.7 & 20.6.8 of
% 'Handbook of Mathematical Functions'
% by Milton Abramowitz and Irene Stegun
% z     - vector of argument values
% m     - integer function order
% q     - scalar parameter value
% ntrms - number of series terms used
% f     - vector of function values

%          by H.W. 10/10/04

if nargin<4, ntrms=50; end
z=z(:); nz=length(z); q2=sqrt(q);
k=0:ntrms-1; sgn=cos(pi*k)';
mh=1+fix(m/2); u1=q2*exp(-z); u2=q2*exp(z);  
if mod(m,2)==0 % even index   
  r=m/2; [a,c]=matue(q,1,1,ntrms); 
  u=c(:,mh); [ubig,s]=max(abs(u));
  p0=(-1)^r/u(s)*sgn;
  if s==1, p0=p0/2; end; s=s-1; 
  f=(bes(k-s,u1).*bes(k+s,u2)+...
     bes(k+s,u1).*bes(k-s,u2))*...
     (p0.*u);
else           % odd index
  [a,c]=matue(q,2,1,ntrms); r1=(m+1)/2;
  u=c(:,mh); [ubig,s]=max(abs(u));
  p0=-(-1)^r1/u(s)*sgn; s=s-1;
  f=(bes(k-s,u1).*bes(k+s+1,u2)+...
     bes(k+s+1,u1).*bes(k-s,u2))*...
     (p0.*u); 
end
end

% function f=bes(n,z)
% % f=bes(n,z) calls function besselj
% % for vector values of n and z.
% % f(i,j)= besselj(n(j),z(i))
% f= besselj(ones(size(z(:)))*n(:)',z(:)*ones(size(n(:)')));