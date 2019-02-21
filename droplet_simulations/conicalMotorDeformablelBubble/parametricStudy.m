%parametric study

%delete(gcp)
%parpool(4)

clear variables
close all

results = '~/Documents/MATLAB/droplet_simulations/results';

%parameters
thetaCone = 2;
nBubble = 80;
Ca = 0.1;
alpha = 0.8;
xStart = 1:4;
dt = 1e-2;
Tend = 500;

parfor l = 1:numel(xStart)

    mainConicalMotorDeformableODEparam(thetaCone,Ca,nBubble,dt,alpha,xStart(l),Tend,results);

end