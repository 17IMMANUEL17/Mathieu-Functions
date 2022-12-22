function v=Ms4v(z,m,q,ntrms)
%  v=Ms4v(z,m,q,ntrms)
if nargin<4, ntrms=50; end
v=Ms1v(z,m,q,ntrms)-1i*Ms2v(z,m,q,ntrms);