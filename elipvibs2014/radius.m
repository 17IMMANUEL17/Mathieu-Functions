function f=radius(a,mi)
%approximation of the peremeter of an ellipse using Ramanujan
%approximation
a_x = a * cosh(mi);
b_x = a * sinh(mi);
tmpH = (a_x-b_x).^2/(a_x+b_x).^2;
perim = pi*(a_x+b_x).*(1+3*tmpH./(10+sqrt(4-3*tmpH)));
f=perim/2/pi;
end