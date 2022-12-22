function [f,args]=cemv(z,m,q,ntrms)
% [f,args]=cemv(z,m,q,ntrms)
% vector form of the modified Mathieu
% function Ce(z,m,q). See 20.6.3 and
% 20.6.4 of Abramowitz&Stegun.

% z     - vector of z values
% m     - vector of function orders
% q     - vector of q values
% ntrms - number of series terms
% f       array in which f(i,j,k)
%         contains cem(z,m,q) eval-
%         uated at z(i),m(j),q(k)
% args    a structure returning the
%         input data in elements
%         args.z, args.m, args.q
if nargin<4, ntrms=100; end
z=z(:); nz=length(z);
m=m(:)'; M=length(m); nq=length(q);
f=zeros(nz,M,nq);
for j=1:M, for k=1:nq
  f(:,j,k)=cem(z,m(j),q(k),ntrms);  
end, end
if nargout==2
  args=struct('z',{z},'m',{m},'q',{q});
end