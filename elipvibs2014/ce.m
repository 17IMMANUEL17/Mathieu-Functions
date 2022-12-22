function v=ce(x,m,q,ntrms)   
% v=ce(x,m,q,ntrms) 
% This function computes the even valued
% Mathieu function ce(x,m,q). See 20.2.3
% of Abramowitz & Stegun.

if nargin<4, ntrms=50; end
if mod(m,2)==0, order=1; n=2*(0:ntrms-1); % even order
else order=2; n=1+2*(0:ntrms-1); end     % odd order 
[a,c]=matue(q,order,1,ntrms);
v=cos(x(:)*n)*c(:,1+fix(m/2)); 