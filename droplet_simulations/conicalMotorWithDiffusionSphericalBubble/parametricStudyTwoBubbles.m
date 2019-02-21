%parametric study

clear variables
close all

results = '~/Documents/MATLAB/droplet_simulations/results';

%parameters
dt = 1e-3;
tol = 1e-3;
ODE = 2;
thetaCone = linspace(pi/64,pi/16,10);
repON = 0.1;
BN = 16;
bubbleCycles = 15;

for i = 1:numel(thetaCone)
    
    display([num2str(i) ' of ' num2str(numel(thetaCone))])

    main_ConicalMotorCambridgeTwoBubblesParamODE(dt,tol,ODE,thetaCone(i),BN,bubbleCycles,results);

end

%exit