%large order asymptotic expansion of J
function f=asy_besy(n,z)
% f=bes(n,z) calls function besselj
% for vector values of n and z.
% f(i,j)= besselj(n(j),z(i))
f=zeros(length(n),length(z));
for i=1:length(n)
    if(n(i)<140)
        f(i,:)= bessely(n(i),z(:));
    else
        f(i,:)=-sqrt(2/(pi*n(i))).*((exp(1)*z)/(2*n(i))).^(-n(i));
    end
end