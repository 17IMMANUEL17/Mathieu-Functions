function figure14_7
% Creates figure 14.7 from 'Computation
% of Special Functions'
x=linspace(0,4,1000); q=10; m=1:4; 
f=squeeze(Ms2v(x,m,q));
plot(x,f(:,1),x,f(:,2),x,f(:,3),...
     x,f(:,4)), grid on
 xlabel('argument x')
 ylabel('Modified Mathieu function Ms2(x,m,q)')
 title('Ms2(x,m,q)  for  q = 10')
 legend('m=1','m=2','m=3','m=4')