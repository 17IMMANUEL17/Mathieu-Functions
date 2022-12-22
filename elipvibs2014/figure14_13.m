function figure14_13
% Creates figure 14.3 from 'Computation
% of Special Functions'
d=0:2:90; x=linspace(0,4,300); q=5; n=1:4; 
f=squeeze(sev(x,n,q));
plot(x,f(:,1),x,f(:,2),x,f(:,3),...
     x,f(:,4)), grid on
 ylabel('Mathieu function Se(x,m,q)')
 title('Se(x,m,q)  for  q = 5')
 legend('Se1','Se2','Se3','Se4')