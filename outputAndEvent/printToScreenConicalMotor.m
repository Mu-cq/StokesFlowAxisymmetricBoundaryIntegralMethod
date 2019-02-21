%print to screen information when starting the simulation

function printToScreenConicalMotor(theta,R,L,x0,alpha,Ca,massFlux)

%output hydrodynamics
display('HYDRODYNAMIC PARAMETERS')
display(['Capillary number is Ca=' num2str(Ca)])
if massFlux==0
    display('Impose volume flux')
elseif massFlux==1
    display('Impose mass flux')
end

%ouput geometry
display('CONE GEOMETRY')
display(['theta=' num2str(theta) ' L=' num2str(L) ' R=' num2str(R)])

%ouput bubble
display('BUBBLE GEOMETRY')
display(['Initial bubble position is x0=' num2str(x0) ' with radius r0=' num2str(alpha)])