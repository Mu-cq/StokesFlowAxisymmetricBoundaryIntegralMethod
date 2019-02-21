%print to screen information when starting the simulation

function printToScreenConicalMotorWithSphere(theta,R,L,x0,alpha,Hcc,beta,endCycle)

%output hydrodynamics
disp('CHEMICAL PARAMETERS')
display(['Hcc=' num2str(Hcc) ' beta=' num2str(beta)])

%ouput geometry
disp('CONE GEOMETRY')
display(['theta=' num2str(theta) ' L=' num2str(L) ' R=' num2str(R)])

%ouput bubble
disp('BUBBLE GEOMETRY')
display(['Initial bubble position is x0=' num2str(x0) ' with radius r0=' num2str(alpha)])

%cycle type
if endCycle==1
    disp('Cycle Type: cycle closes when CM bubble exits the cone')
elseif endCycle==2
    disp('Cycle Type: bubble is eliminated when shrinking')
elseif endCycle==3
    disp('Cycle Type: largest bubble is eliminated when third bubble nucleates')
end