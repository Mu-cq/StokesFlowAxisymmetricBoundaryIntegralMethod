%parametric study

%delete(gcp)
%parpool(4)

clear variables
close all

%results = '~/Documents/MATLAB/droplet_simulations/server';
results = '~/Documents/MATLAB/droplet_simulations/results/forThesis';

%parameters
nDrop = round(logspace(1,2.3,10));
Ca = 0.1;
Tend = 100;
dt =  [0.1*ones(1,6) 0.02*ones(1,1) 0.01*ones(1,3)];
ODE = 2;
lambda = 0.5;
curved = 0;

if curved==1
    nameShort = 'SPLINE';
elseif curved==0
    nameShort = 'STRAIGTH';
else
    nameShort = '';
end

for l = 1:numel(nDrop)
    
    display(['Simulation ' num2str(l) ' of ' num2str(numel(nDrop))])

    main_oneDropBEMparam(nDrop(l),Ca,lambda,dt(l),Tend,curved,results,nameShort);

end

%exit