function v=Mc4v(z,m,q,ntrms)
%  v=Mc4v(z,m,q,ntrms)
if nargin<4, ntrms=50; end
v=Mc1v(z,m,q,ntrms)-1i*Mc2v(z,m,q,ntrms);