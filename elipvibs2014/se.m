function v=se(x,m,q,ntrms)   
% v=se(x,m,q,ntrms) 
% This function computes the odd valued
% Mathieu function se(x,m,q). 
% See 20.2.4 of Abramowitz & Stegun.

if m<=0
 error('In se(x,m,q), order m <= 0 is not allowed.')
end
if nargin<4, ntrms=50; end
if mod(m,2)==0
  order=1; n=2*(1:ntrms);     % even order
  [a,c]=matue(q,order,2,ntrms); 
  v=sin(x(:)*n)*c(:,m/2);
else
  order=2; n=1+2*(0:ntrms-1); % odd order
  [a,c]=matue(q,order,2,ntrms);
  v=sin(x(:)*n)*c(:,1+fix(m/2));
end
end