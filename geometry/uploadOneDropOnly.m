%upload initial bubble shape

function [x,y,PARAMnew] = uploadOneDropOnly(ODE,Tend,maxDT,PARAMnew)

filename = ['oneDropBEM_ODE=' num2str(ODE) '_n=' num2str(PARAMnew.nUP) '_BC=' num2str(PARAMnew.typeBCstokes) '_Ca=' num2str(PARAMnew.CaUP) '_visc=' num2str(PARAMnew.lambdaUP) '_D=' num2str(PARAMnew.D) '_maxDT=' num2str(maxDT) '_Tend=' num2str(Tend) '.mat'];
cd(PARAMnew.res)
load(filename)
cd(PARAMnew.here)

%select last shape
Y = allRes{2};
Ylast = Y{end};
x = Ylast(1:2:end-1);
y = Ylast(2:2:end);

PARAMnew.n(1) = numel(x)-1;