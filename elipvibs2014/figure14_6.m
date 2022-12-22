function figure14_6
% Creates figure 14.6 from 'Computation
% of Special Functions'
x=linspace(0,4,1000); q=10; m=0:3; 
f=squeeze(Mc2v(x,m,q));
plot(x,f(:,1),x,f(:,2),x,f(:,3),...
     x,f(:,4)), grid on
 xlabel('argument x')
 ylabel('Modified Mathieu function Mc2(x,m,q)')
 title('Mc2(x,m,q)  for  q = 10')
 legend('m=0','m=1','m=2','m=3')