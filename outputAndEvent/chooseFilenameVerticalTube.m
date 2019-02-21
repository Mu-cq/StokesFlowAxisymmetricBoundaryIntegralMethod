%choose filename

function filename = chooseFilenameVerticalTube(PARAM,nStart,L,Tend,dt,alpha,x0,Bond,ODE,lambda,rep)

%check if saving directory exist
here = pwd;
cd(PARAM.res)
cd(here)

%filename = ['verticalTube_ODE=' num2str(ODE) '_x0=' num2str(x0) '_alpha=' num2str(alpha) '_Bond=' num2str(Bond) '_lambda=' num2str(lambda) '_L=' num2str(L) '_nStart=' num2str(nStart) '_Tend=' num2str(Tend) '_dt=' num2str(dt) '.mat'];
filename = ['verticalTube_ODE=' num2str(ODE) '_rep=' num2str(rep) '_x0=' num2str(x0) '_alpha=' num2str(alpha) '_Bond=' num2str(Bond) '_lambda=' num2str(lambda) '_L=' num2str(L) '_nStart=' num2str(nStart) '_Tend=' num2str(Tend) '_dt=' num2str(dt) '.mat'];