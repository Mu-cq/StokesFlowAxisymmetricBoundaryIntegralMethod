%print to screen information when starting the simulation

function printToScreenBubbleCloseToWall(gap,Ca,lambda)

%output hydrodynamics
disp('HYDRODYNAMIC PARAMETERS')
disp(['Capillary number is Ca=' num2str(Ca) ', Viscosity ratio is lambda=' num2str(lambda)])

%ouput geometry
disp('GEOMETRY')
disp(['Initial gap size is gap=' num2str(gap)])