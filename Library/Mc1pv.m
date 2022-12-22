function [fp,args]=Mc1pv(z,m,q,ntrms)
% [fp,args]=Mc1pv(z,m,q,ntrms)

% 3D Array form of the derivative with
% respect to z of the modified Mathieu
% function Mc1(z,m,q) defined by 
% 20.6.9 & 20.6.10 of
% 'Handbook of Mathematical Functions'
% by Milton Abramowitz and Irene Stegun
% z     - vector of argument values
% m     - integer function order
% q     - scalar parameter value
% ntrms - number of series terms used
% fp    - vector of function values
% z     - vector of z values
% m     - vector of function orders
% q     - vector of q values
% ntrms - number of series terms
% fp      an array in which fp(i,j,k) is
%         the derivative of Ms1(z,m,q)
%         evaluated at z(i),m(j),q(k)
% args    a structure returning the
%         input data in elements
%         args.z, args.m, args.q

%         by H.W. 10/10/04
if nargin<4, ntrms=50; end
z=z(:); nz=length(z); m=m(:)';
M=length(m); nq=length(q); 
fp=zeros(nz,M,nq);
for j=1:M, for k=1:nq
  fp(:,j,k)=Mc1p(z,m(j),q(k),ntrms);  
end, end
if nargout==2
  args=struct('z',{z},'m',{m},'q',{q});
end