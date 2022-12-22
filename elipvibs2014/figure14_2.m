function figure14_2
% Creates figure 14.2 from 'Computation
% of Special Functions'
d=0:2:90; x=linspace(0,10,300); q=10; n=0:3; 
f=squeeze(cev(x,n,q));
plot(x,f(:,1),x,f(:,2),x,f(:,3),...
     x,f(:,4)), grid on
 
 ylabel('Mathieu function ce(x,m,q)')
 title('ce(x,m,q)  for  q = 10')
 legend('ce0','ce1','ce2','ce3')