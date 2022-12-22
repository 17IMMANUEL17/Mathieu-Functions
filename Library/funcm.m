
function v=funcm(w,z,n,h2,ftype,ntrms)
% v=funcm(w,z,n,h2,ftype,ntrms)
if nargin<6, ntrms=50; end; q=(h2*w)^2; 
switch ftype
  case 1, v=Mc1(z,n,q,ntrms);
  case 2, v=Ms1(z,n,q,ntrms);
  case 3, v=Mc1p(z,n,q,ntrms); %disp([w,z,n,ftype]),pause
  case 4, v=Ms1p(z,n,q,ntrms);
end
v=sign(v).*(abs(v)).^(1/5);