function figure14_8
% Creates figure 14.8 from 'Computation
% of Special Functions'
x=1.317; q=.01:.1:200; m=0:3; 
f=squeeze(Mc1v(x,m,q));
plot(q,f(1,:),q,f(2,:),q,f(3,:),...
     q,f(4,:)), grid on
 xlabel('argument q')
 ylabel('Modified Mathieu function Mc1(x,m,q)')
 title('Mc1(x,m,q)  for  x = 1.317')
 legend('m=0','m=1','m=2','m=3')