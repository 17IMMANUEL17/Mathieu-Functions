function figure14_10
% Creates figure 14.10 from 'Computation
% of Special Functions'
x=1.317; q=.1:.1:20; m=0:3; 
f=squeeze(Mc2v(x,m,q));
plot(q,f(1,:),q,f(2,:),q,f(3,:),...
     q,f(4,:)), grid on
 axis([0,20,-.6,.6])
 xlabel('argument q')
 ylabel('Modified Mathieu function Mc2(x,m,q)')
 title('Mc2(x,m,q)  for  x = 1.317')
 legend('m=0','m=1','m=2','m=3')