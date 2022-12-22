function f=ldaC_S(n,mu,q,nMax)
% depending on q and mu we can get to a larger order without any
% instability and so we have to use the asymptotic approximation presented
% in the last chapter of the thesis after a precise n (that depends, as
% just said, on q and mu). This is the approximation for the C eigenvalues of the single layer.
if(n <= 1700)
    f= 1i*pi/2 * Mc1(mu,n,q,nMax) *( Mc1(mu,n,q,nMax) + 1i * Mc2(mu,n,q,nMax));
else
    f= 1i*pi/2*1/n*0.317836;
end

end