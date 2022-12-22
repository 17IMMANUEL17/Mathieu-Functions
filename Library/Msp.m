function fp=Msp(z,m,q,type,ntrms)
% fp=Msp(z,m,q,type,ntrms)

% Derivative with respect to z of the
% modified Mathieu function Ms(z,m,q,type)
% defined by 20.6.9 & 20.6.10 of
% 'Handbook of Mathematical Functions'
% by Milton Abramowitz and Irene Stegun
% z     - vector of argument values
% m     - integer function order
% q     - scalar parameter value
% type  - 1,2,3, or 4 to use Besselj,
%         bessely,hankel1, or hankel2
% ntrms - number of series terms used
% fp    - vector of function values

%       by H.W. 10/10/04

if nargin<5, ntrms=50; end
if nargin<4, type=1; end
z=z(:); z1=sqrt(q)*exp(-z); 
z2=sqrt(q)*exp(z); w=ones(1,ntrms);
z1p=-z1*w; z2p=z2*w;
if mod(m,2)==0 % even order
  k=1:ntrms; [a,c]=matue(q,1,2,ntrms);
  u=c(:,m/2); r=m/2;
  [ubig,s]=max(abs(u)); p0=(-1)^r/u(s);
  %f=(bes(k-s,z1).*bes(k+s,z2)-...
  %   bes(k+s,z1).*bes(k-s,z2))*...
  %   (p0*cos(k*pi)'.*u);
  
  %f1p=z1p.*(besp(k-s,z1).*bes(k+s,z2)-...
  %          besp(k+s,z1).*bes(k-s,z2));
  %f2p=z2p.*(bes(k-s,z1).*besp(k+s,z2)-...
  %          bes(k+s,z1).*besp(k-s,z2));
    
  f1p=z1p.*(besp(k-s,z1).*Z(k+s,type,z2)-...
            besp(k+s,z1).*Z(k-s,type,z2));
  f2p=z2p.*(bes(k-s,z1).*Zp(k+s,type,z2)-...
            bes(k+s,z1).*Zp(k-s,type,z2));
  fp=(f1p+f2p)*(p0*cos(k*pi)'.*u);     
  
else           % odd order
  k=0:ntrms-1; [a,c]=matue(q,2,2,ntrms); 
  u=c(:,1+fix(m/2)); [ubig,s]=max(abs(u));
  
  r=(m-1)/2; p0=(-1)^r/u(s); s=s-1;
  
  %r=(m-1)/2; p0=(-1)^r/u(s);
  
  %f=(bes(k-s,z1).*bes(k+s+1,z2)-...
  %   bes(k+s+1,z1).*bes(k-s,z2))*...
  %   (p0*cos(k*pi)'.*u);
  
  %f1p=z1p.*(besp(k-s,z1).*bes(k+s+1,z2)-...
  %          besp(k+s+1,z1).*bes(k-s,z2));
  %f2p=z2p.*(bes(k-s,z1).*besp(k+s+1,z2)-...
  %          bes(k+s+1,z1).*besp(k-s,z2));
        
  f1p=z1p.*(besp(k-s,z1).*Z(k+s+1,type,z2)-...
            besp(k+s+1,z1).*Z(k-s,type,z2));
  f2p=z2p.*(bes(k-s,z1).*Zp(k+s+1,type,z2)-...
            bes(k+s+1,z1).*Zp(k-s,type,z2));          
  fp=(f1p+f2p)*(p0*cos(k*pi)'.*u);      
end

%========================================

% function f=bes(n,z)
% % f=bes(n,z) computes besselj for vector
% % values of n and z
% f= besselj(ones(size(z(:)))*n(:)',z(:)*ones(size(n(:)')));
% 
% %========================================
% 
% function fp=besp(n,z)
% % fp=besp(n,z) computes the derivative of 
% % besselj for vector values of n and z
% n=n(:)'; z=z(:);
% fp=(bes(n-1,z)-bes(n+1,z))/2;