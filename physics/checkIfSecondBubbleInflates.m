%Choose bubble position and size

function status = checkIfSecondBubbleInflates(t,var,fCompute2bubbles,z0,maxPHI,beta,Hcc)

%comoute bubble radius
r0 = 2*Hcc/(maxPHI-beta*Hcc);
r0 = r0+0.2*r0;

%compute flow rates
[~,~,~,QmassRight,QmassLeft,~] = fCompute2bubbles(t,[var(1); var(2); z0; var(3); r0]);

%if second bubble (left) inflates, stop simulation
% if QmassLeft>0
%     
%     display('New bubble nucleates and grows')
%     status = 1;
%     
% else
%     
%     display('New bubble nucleates but it does not grow')
%     status = 0;
%     
% end
if QmassRight<=0
    
    display('New bubble nucleates and the old one stop growing')
    status = 1;
    
else
    
    display('New bubble nucleates but the large one still grows')
    status = 0;
    
end









