function figure14_11
% Creates figure 14.11 from 'Computation
% of Special Functions'
x=1.317; q=.1:.1:20; m=1:4; 
f=squeeze(Ms2v(x,m,q));
plot(q,f(1,:),q,f(2,:),q,f(3,:),...
     q,f(4,:)), grid on
 axis([0,20,-.6,.6])
 xlabel('argument q')
 ylabel('Modified Mathieu function Ms2(x,m,q)')
 title('Ms2(x,m,q)  for  x = 1.317')
 legend('m=1','m=2','m=3','m=4')