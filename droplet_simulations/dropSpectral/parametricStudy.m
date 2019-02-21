%parametrci study

delete(gcp)
parpool(3)

clear variables
close all

%results folder
res = '~/Documents/MATLAB/droplet_simulations/results';

%parameters
visc = 1;
Ca = 0;
elem = 50;
Thorizon = 100:109;
dt = 5e-2;
A0 = [];
upShape = 4;
legendre = 0;
%Ddef = [0.007 0.0162 0.0377 0.0874 0.2026];
Ddef = 0;
%CaUP = [0.0046 0.0105 0.024 0.055 0.115 0.21];
%CaUP = [0.004 0.009 0.02 0.044 0.081];
CaUP = 0.125;

%parameters upload
ThorizonUp = [];
A0up = [];

parfor l = 1:numel(Thorizon)

    %run simulation
    mainDropSpectralParametric(A0,Thorizon(l),Ca,visc,elem,ThorizonUp,A0up,dt,Ddef,CaUP,upShape,legendre,res);

end

%exit