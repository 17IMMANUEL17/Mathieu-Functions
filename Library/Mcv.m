function [f,args]=Mcv(z,m,q,type,ntrms)
% [f,args]=Mcv(z,m,q,type,ntrms)

% 3D array form of the modified Mathieu
% function Mc(z,m,q) defined by 
% 20.4.7, 20.6.7 & 20.6.8 of
% 'Handbook of Mathematical Functions'
% by Milton Abramowitz and Irene Stegun
% z     - vector of argument values
% m     - integer function order
% q     - scalar parameter value
% type  - 1,2,3, or 4 to use Besselj,
%         bessely,hankel1, or hankel2
% ntrms - number of series terms used
% f       array in which f(i,j,k)
%         contains Mc1(z,m,q) eval-
%         uated at z(i),m(j),q(k)
% args    a structure returning the
%         input data in elements
%         args.z, args.m, args.q

%          by H.W. 10/10/04

if nargin<5, ntrms=50; end
if nargin<4, type=1; end
z=z(:); nz=length(z);
m=m(:)'; M=length(m); nq=length(q);
f=zeros(nz,M,nq);
for j=1:M, for k=1:nq
  f(:,j,k)=Mc(z,m(j),q(k),type,ntrms);  
end, end
if nargout==2
  args=struct(...
      'z',{z},'m',{m},'q',{q},'type',{type});
end