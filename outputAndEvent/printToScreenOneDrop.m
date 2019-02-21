%print to screen information when starting the simulation

function printToScreenOneDrop(Ca,Bond,lambda,dropFrame,algo,loopEdge,totalLoopEdge)

%algorithm
disp('ALGORITHM')
if algo==1
    disp('DNS simulation')
elseif algo==2
    disp('Minimal seed optimization')
elseif algo==3
    disp(['Edge tracking, loop number ' num2str(loopEdge) ' of ' num2str(totalLoopEdge)])
else
   error('Not implemented') 
end

%output hydrodynamics
disp('HYDRODYNAMIC PARAMETERS')
if isempty(Bond)==1 && isempty(Ca)==0
    display(['Capillary number is Ca=' num2str(Ca) ', Viscosity ratio is lambda=' num2str(lambda)])
elseif isempty(Bond)==0 && isempty(Ca)==1
    display(['Bond number is Bo=' num2str(Bond) ', Viscosity ratio is lambda=' num2str(lambda)])
else
    error('You have to choose Bond Or Capillary number')
end

%frame of reference
if dropFrame==0
    frame = 'lab frame';
elseif dropFrame==1
    frame = 'drop frame';
else
    error('Not implemented')
end
disp(['Simulation is in the ' frame])