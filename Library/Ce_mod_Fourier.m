function v=Ce_mod_Fourier(x,m,q,ntrms)   
% v=Ce(x,m,q,ntrms) Ce obtained through the Fourier expansion

if nargin<4, ntrms=50; end
if mod(m,2)==0, order=1; n=2*(0:ntrms-1); % even order
else order=2; n=1+2*(0:ntrms-1); end     % odd order 
[a,c]=matue(q,order,1,ntrms);
v=cosh(x(:)*n)*c(:,1+fix(m/2)); 