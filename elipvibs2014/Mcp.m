function fp=Mcp(z,m,q,type,ntrms)
% fp=Mcp(z,m,q,type,ntrms)

% Derivative with respect to z of the
% modified Mathieu function Mc(z,m,q,type)
% defined by 20.4.7, 20.6.7 & 20.6.8 of
% 'Handbook of Mathematical Functions'
% by Milton Abramowitz and Irene Stegun
% z     - vector of argument values
% m     - integer function order
% q     - scalar parameter value
% type  - 1,2,3,4 for modified Mathieu
%         functions of kinds 1,2,3, or 4
% ntrms - number of series terms used
% fp    - vector of derivative values

%          by H.W. 3/20/06

if nargin<5, ntrms=50; end
if nargin<4, type=1; end
z=z(:); nz=length(z); q2=sqrt(q);
k=0:ntrms-1; sgn=cos(pi*k)'; w=ones(1,ntrms);
mh=1+fix(m/2); u1=q2*exp(-z); u2=q2*exp(z); 
u1p=-u1*w; u2p=u2*w; 
if mod(m,2)==0 % even index   
  r=m/2; [a,c]=matue(q,1,1,ntrms); 
  u=c(:,mh); [ubig,s]=max(abs(u));
  p0=(-1)^r/u(s)*sgn;
  if s==1, p0=p0/2; end; s=s-1; 
  %f=(besselj(k-s,u1).*besselj(k+s,u2)+...
  %   besselj(k+s,u1).*besselj(k-s,u2))*...
  %   (p0.*u);
  
  %f1p=u1p.*(besp(k-s,u1).*bes(k+s,u2)+...
  %          besp(k+s,u1).*bes(k-s,u2));
  %f2p=u2p.*(bes(k-s,u1).*besp(k+s,u2)+...
  %          bes(k+s,u1).*besp(k-s,u2));
        
  f1p=u1p.*(besp(k-s,u1).*Z(k+s,type,u2)+...
            besp(k+s,u1).*Z(k-s,type,u2));
  f2p=u2p.*(bes(k-s,u1).*Zp(k+s,type,u2)+...
            bes(k+s,u1).*Zp(k-s,type,u2)); 
  fp=(f1p+f2p)*(p0.*u); 
  
else           % odd index
  [a,c]=matue(q,2,1,ntrms); r1=(m+1)/2;
  u=c(:,mh); [ubig,s]=max(abs(u));
  p0=-(-1)^r1/u(s)*sgn; s=s-1;
  %f=(besselj(k-s,u1).*besselj(k+s+1,u2)+...
  %   besselj(k+s+1,u1).*besselj(k-s,u2))*...
  %   (p0.*u); 
  
  %f1p=u1p.*(besp(k-s,u1).*bes(k+s+1,u2)+...
  %          besp(k+s+1,u1).*bes(k-s,u2));
  %f2p=u2p.*(bes(k-s,u1).*besp(k+s+1,u2)+...
  %          bes(k+s+1,u1).*besp(k-s,u2));
        
  f1p=u1p.*(besp(k-s,u1).*Z(k+s+1,type,u2)+...
            besp(k+s+1,u1).*Z(k-s,type,u2));
  f2p=u2p.*(bes(k-s,u1).*Zp(k+s+1,type,u2)+...
            bes(k+s+1,u1).*Zp(k-s,type,u2));          
  fp=(f1p+f2p)*(p0.*u);      
end

%========================================

function f=bes(n,z)
% f=bes(n,z) calls function besselj
% for vector values of n and z.
% f(i,j)= besselj(n(j),z(i))
f= besselj(ones(size(z(:)))*n(:)',z(:)*ones(size(n(:)')));

%========================================

function fp=besp(n,z)
% fp=besp(n,z) computes the derivative of 
% besselj for vector values of n and z
n=n(:)'; z=z(:);
fp=(bes(n-1,z)-bes(n+1,z))/2;