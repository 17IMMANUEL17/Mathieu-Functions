function v=Ms4pv(z,m,q,ntrms)
%  derivative of Ms4 with respect to z
if nargin<4, ntrms=50; end
v=Ms1pv(z,m,q,ntrms)-i*Ms2pv(z,m,q,ntrms);