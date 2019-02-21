%choose filename

function filename = chooseFilenameConicalMotorCycle(PARAM,nPerLenght,beta,K,BN,L,Tend,thetaCone,dt,ODE,tol,bubbleCycles)

%check if saving directory exist
here = pwd;
cd(PARAM.res)
cd(here)

filename = ['coneWithDiffusion_nCycle=' num2str(bubbleCycles) '_ODE=' num2str(ODE) '_tol=' num2str(tol) '_beta=' num2str(beta) '_Hcc=' num2str(K) '_BN=' num2str(BN) '_L=' num2str(L) '_nPerL=' num2str(nPerLenght) '_rep=' num2str(sum(PARAM.repulsiveForces)) '_coeffRep=' num2str(PARAM.coeffRepulsive/sin(PARAM.rotate(4))) '_distRep=' num2str(PARAM.repulsiveOn) '_theta=' num2str(thetaCone) '_Tend=' num2str(Tend) '_dt=' num2str(dt) '.mat'];