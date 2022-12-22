function [a,c,v]=matue(q,order,ftype,ntrms)
% [a,c,v]=matue(q,order,ftype,ntrms)

% For a given value of q, this function computes
% values of parameter a in the Mathieu equation 
%     y"(x)+(a-2*q*cos(2*x))*y(x)=0
% The coefficients in the Fourier series
% expansions of the Mathieu functions
% ce(x,n,q) and se(x,n,q) are also obtained.
% q     - parameter in the Mathieu equation
%         proportional to the square of the
%         related natural frequency of an
%         elliptic membrane. This value can
%         be complex.
% order - 1 for even index, 2 for odd index
% ftype - 1 for even valued, 2 for odd-valued
% ntrms - number of terms used to approximate
%         the Fourier expansion of ce or se
% a     - eigenvalues which produce periodic
%         solutions of the Mathieu equation
% c     - matrix of eigenvectors defining the
%         Fourier coefficients of the cosine
%         or sine series
% v     - tridiagonal matrix which determines 
%         the eigenvalues and eigenvectors

%         HBW 12/01/04

if nargin<4, ntrms=50; end
if nargin==0, q=10; order=1; ftype=1; end
v=diag(q*ones(ntrms-1,1),1); v=v+v.';
if ftype==1;  % even valued
  if order==1 % even order
    v=v+diag((2*(0:ntrms-1)).^2);
    v(1,2)=sqrt(2)*q; v(2,1)=v(1,2);     
  else        % odd order
    v=v+diag((1+2*(0:ntrms-1)).^2);
    v(1,1)=1+q;
  end
else % ftype==2, odd valued
  if order==1  % even order
    v=v+diag((2*(1:ntrms)).^2);
  else         % odd order 
    v=v+diag((1+2*(0:ntrms-1)).^2);  
    v(1,1)=1-q;
  end
end
if abs(q)>1000,v=v/q; end
[c,a]=eig(v); 
if abs(q)>1000, a=a*q; end
a=diag(a); [~,k]=sort(real(a));
a=a(k); c=c(:,k);
if order==1 && ftype==1
  c(1,:)=c(1,:)/sqrt(2);
end    
% Normalize the series coefficients
cv=sum(conj(c).*c);
if order==1 && ftype==1
  cv=cv+conj(c(1,:)).*c(1,:);
end
if~isreal(q) % Complex eigenvector case
  c=ones(ntrms,1)*(1./sqrt(cv)).*c;
  return
end
% Case where q is real. Make ce(0,q)
% or se'(0,q) positive.

c=c.*cv(ones(ntrms,1),:); 