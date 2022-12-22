function figure14_12
d=0:2:90; x=linspace(0,3,300); q=1; n=0:3; 
f=squeeze(Cev_mod_Fourier(x,n,q));
plot(x,f(:,1),x,f(:,2),x,f(:,3),...
     x,f(:,4)), grid on
 
 ylabel('Mathieu function Ce(x,m,q)')
 title('Ce(x,m,q)  for  q = 1')
 legend('Ce0','Ce1','Ce2','Ce3')