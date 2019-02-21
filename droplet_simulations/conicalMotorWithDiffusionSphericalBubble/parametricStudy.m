%parametric study

%parpool(16)

clear variables
close all

results = '~/Documents/MATLAB/droplet_simulations/results';

%parameters
dt = 1e-3;
tol = 1e-3;
ODE = 2;
thetaCone = linspace(pi/64,pi/6,10);
%thetaCone = thetaCone([8 4 1]);
thetaCone = 6/180*pi;
%L = 10;
L = 10;
repON = 0.1;
repIntensity = 1e7;
%BN = [16:19  21:25];
BN = [];
Hcc = 10;
beta = 10;
cycleType = 1;
Tend = 1;
checkpoint = 1e3;

for l = 1:numel(Hcc)
    for k = 1:numel(beta)
        for i = 1:numel(thetaCone)

            display([num2str(i) ' of ' num2str(numel(thetaCone))])

            main_ConicalMotorCambridgeODEparam(dt,tol,ODE,thetaCone(i),L,repON,repIntensity,BN,Hcc(l),beta(k),cycleType,Tend,checkpoint,results);

        end
    end
end

%exit