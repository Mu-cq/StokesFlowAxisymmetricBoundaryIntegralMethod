%choose filename

function filename = chooseFilenameBubbleCloseToWall(PARAM,nStart,Tend,dt,alpha,Ca,gap,ODE,lambda)

%check if saving directory exist
here = pwd;
cd(PARAM.res)
cd(here)

filename = ['bubbleCloseToWall_ODE=' num2str(ODE) '_gap=' num2str(gap) '_alpha=' num2str(alpha) '_Ca=' num2str(Ca) '_lambda=' num2str(lambda) '_nStart=' num2str(nStart) '_Tend=' num2str(Tend) '_dt=' num2str(dt) '.mat'];