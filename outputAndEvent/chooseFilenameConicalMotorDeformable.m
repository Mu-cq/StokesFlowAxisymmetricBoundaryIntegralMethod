%choose filename

function filename = chooseFilenameConicalMotorDeformable(PARAM,nStart,L,Tend,thetaCone,dt,alpha,x0,Ca,ODE,massFlux)

%check if saving directory exist
here = pwd;
cd(PARAM.res)
cd(here)

filename = ['coneBubbleDeformable_massFlux=' num2str(massFlux) '_ODE=' num2str(ODE) '_x0=' num2str(x0) '_alpha=' num2str(alpha) '_Ca=' num2str(Ca) '_L=' num2str(L) '_nStart=' num2str(nStart) '_rep=' num2str(sum(PARAM.repulsiveForces)) '_coeffRep=' num2str(PARAM.coeffRepulsive) '_distRep=' num2str(PARAM.repulsiveOn) '_theta=' num2str(thetaCone) '_Tend=' num2str(Tend) '_dt=' num2str(dt) '.mat'];