%parametric study

%delete(gcp)
%parpool(4)

clear variables
close all

%results = '~/Documents/MATLAB/droplet_simulations/server';
results = '~/Documents/MATLAB/droplet_simulations/results';

%parameters
nDrop = 20;
nUpload = nDrop+110;
Bond = flip(0.2:0.05:0.9);
BondUpload = [1 Bond(1:end-1)];
uploadRes = 0;
alpha = 1.2;
Tend = 1000;
dt =  1;
ODE = 2;
lambda = 1;
ODEtol = 1e-3;
volTol = 1e-3;
dtUP = dt;
remeshTypeBubble = 4;
adaptLevel = 1.2;

for l = 1:numel(Bond)
    
    display(['Simulation ' num2str(l) ' of ' num2str(numel(Bond))])

    main_verticalPipeDropParam(alpha,nDrop,nUpload,Bond(l),lambda,BondUpload(l),results,uploadRes,Tend,dt,dtUP,ODE,volTol,ODEtol,adaptLevel,remeshTypeBubble);

end

%exit