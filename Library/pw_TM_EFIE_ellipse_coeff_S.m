function pw_coeff = pw_TM_EFIE_ellipse_coeff_S(phi_i, n,mu,q,nMax)
if n == 0
    pw_coeff =0;
else if n < 0
    pw_coeff = (1i).^(n).*Ms1(mu,abs(n),q,nMax).*exp(1i.*n.*phi_i);
else
    pw_coeff = (1i).^(n).*Ms1(mu,n,q,nMax).*exp(1i.*n.*phi_i);
end
end
