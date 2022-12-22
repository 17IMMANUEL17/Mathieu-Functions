function figure14_9
% Creates figure 14.9 from 'Computation
% of Special Functions'
x=1.317; q=0:.1:20; m=1:4; 
f=squeeze(Ms1v(x,m,q));
plot(q,f(1,:),q,f(2,:),q,f(3,:),...
     q,f(4,:)), grid on
 xlabel('argument q')
 ylabel('Modified Mathieu function Ms1(x,m,q)')
 title('Ms1(x,m,q)  for  x = 1.317')
 legend('m=1','m=2','m=3','m=4')