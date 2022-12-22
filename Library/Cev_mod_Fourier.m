function [f,args]=Cev_mod_Fourier(z,m,q,ntrms)
% [f,args]=cev(z,m,q,ntrms)
% vector form of function Ce(z,m,q)

% z     - vector of z values
% m     - vector of function orders
% q     - vector of q values
% ntrms - number of series terms
% f       array in which f(i,j,k)
%         contains ce(z,m,q) eval-
%         uated at z(i),m(j),q(k)
% args    a structure returning the
%         input data in elements
%         args.z, args.m, args.q
if nargin<4, ntrms=100; end
z=z(:); nz=length(z);
m=m(:)'; M=length(m); nq=length(q);
f=zeros(nz,M,nq);
for j=1:M, for k=1:nq
  f(:,j,k)=Ce_mod_Fourier(z,m(j),q(k),ntrms);  
end, end
if nargout==2
  args=struct('z',{z},'m',{m},'q',{q});
end