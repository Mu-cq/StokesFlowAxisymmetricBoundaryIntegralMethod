%upload initial bubble shape

function [x,y,PARAMnew] = uploadBubbleShape(Bo,Ode,x0,alpha,lambda,L,nOldDrop,Tend,dt,rep,PARAMnew)

disp(['Initial condition is from previous simulation with Bo=' num2str(Bo) ' and lambda=' num2str(lambda)])

filename = ['verticalTube_ODE=' num2str(Ode) '_rep=' num2str(rep) '_x0=' num2str(x0) '_alpha=' num2str(alpha) '_Bond=' num2str(Bo) '_lambda=' num2str(lambda) '_L=' num2str(L) '_nStart=' num2str(nOldDrop) '_Tend=' num2str(Tend) '_dt=' num2str(dt) '.mat'];
cd(PARAMnew.res)
load(filename)
cd(PARAMnew.here)

%select last shape
Ylast = Y{end};
x = Ylast(1:2:end-1);
y = Ylast(2:2:end);

PARAMnew.n(4) = numel(x)-1;