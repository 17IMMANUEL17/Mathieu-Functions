function v=Mc4pv(z,m,q,ntrms)
%derivative of Mc4 with respect to z
if nargin<4, ntrms=50; end
v=Mc1pv(z,m,q,ntrms)-1i*Mc2pv(z,m,q,ntrms);