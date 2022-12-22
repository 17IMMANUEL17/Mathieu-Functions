function [q,a,b]= figure14_1
% Creates Figure 14.1 of 'Computation
% of Special Functions'
q=(-30:.25:30); nmax=6; ntrms=50;
a=eigmatu(q,nmax,1,ntrms);
b=eigmatu(q,nmax,2,ntrms);
r=[-30,30,-12,42]; 
plot(q,a,'k-'), axis(r),hold on 
plot(q,b,'k:'), axis(r)
xlabel('q axis')
ylabel('eigenvalues of ce and se')
title('Eigenvalues of Mathieu Functions')
xg=[684 1264 518 877 1513 1970 1264 ...
    -1887 -2205 -1264 -504 -421 -311]/100;
yg=[-607 -339 576 1255 2344 3465 4002 ...
    3876 3576 1855 1318 402 118]/100;
lab=['a0';'a1';'a2';'a3';'a4';'a5';...
      'a6';'b6';'b5';'b4';'b3';'b2';'b1'];
for j=1:13
  text(xg(j),yg(j),lab(j,:));
end   
hold off, shg

%==============================

function f=eigmatu(q,nmax,ftype,ntrms)
% f=eigmatu(q,nmax,ftype,ntrms)

% This function returns an array of
% characteristic values for the 
% angular Mathieu functions ce or se.
% q     - a vector of q values
% nmax  - highest function order used.
%         Columns of f correspond to
%         0:nmax for ce or 1:nmax for 
%         se.
% ftype - 1 for ce, 2 for se
% ntrms - number of series terms used
% f     - a matrix of eigenvalues with
%         q varying by row index and
%         order vary by column index.
if nargin<4, ntrms=50; end
if ftype==1, f=eigc(q,nmax,ntrms);
else, f=eigs(q,nmax,ntrms); end

%==============================

function f=eigc(q,nmax,ntrms)
% f=eigc(q,mmax,ntrms)
% characteristic values for ce
if nargin<3, ntrms=50; end
nq=length(q); nh=1:1+ceil(nmax/2);
f=zeros(nq,nmax+1); ncol=1:nmax+1;
for j=1:nq
  ae=matue(q(j),1,1,ntrms);
  ao=matue(q(j),2,1,ntrms);
  t=[ae(nh),ao(nh)]'; 
  t=t(:)'; f(j,:)=t(ncol);
end

%==============================

function f=eigs(q,nmax,ntrms)
% f=eigs(q,mmax,ntrms)
% characteristic values for se
if nargin<3, ntrms=50; end
nq=length(q); nh=1:ceil(nmax/2);
f=zeros(nq,nmax); ncol=1:nmax;
for j=1:nq
  be=matue(q(j),1,2,ntrms);
  bo=matue(q(j),2,2,ntrms);
  t=[bo(nh),be(nh)]'; 
  t=t(:)'; f(j,:)=t(ncol);
end  