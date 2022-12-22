function figure14_3
% Creates figure 14.3 from 'Computation
% of Special Functions'
d=0:2:90; x=linspace(0,10,300); q=10; n=1:4; 
f=squeeze(sev(x,n,q));
plot(x,f(:,1),x,f(:,2),x,f(:,3),...
     x,f(:,4)), grid on
 ylabel('Mathieu function se(x,m,q)')
 title('se(x,m,q)  for  q = 10')
 legend('se1','se2','se3','se4')