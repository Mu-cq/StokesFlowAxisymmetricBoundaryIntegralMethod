%activate spectral remesh

function [z0,r0] = findLocationAndRadiusOfSecondBubble(t,var,fCompute,beta,Hcc,PARAM)
    
%compute concentration along the axis
[~,yLaplace,~,xVel,yVel,PARAMvel] = fCompute(t,var);
Xline = linspace(-10+var(1),20+var(1),200);
Yline = zeros(1,numel(Xline));
[~,~,PHIfield] = computeConcentrationField(Xline,Yline,xVel,yVel,yLaplace,PARAMvel,0,1);

%find position and radius new bubble
[maxPHI,indPHI] = max(PHIfield);
z0 = Xline(indPHI);
r0 = 2*Hcc/(maxPHI-beta*Hcc);
r0 = r0+0.1;

%check if bubble radius is too large
if r0>PARAM.maxRadiusNucleation
    error(['Initial radius it too large, r=' num2str(r0)])
end

if r0<PARAM.smallestRadius
    r0 = PARAM.smallestRadius;
    display(['Minimum bubble size allowed is r=' num2str(PARAM.smallestRadius)])
end

%check if bubble radius is negative
if r0==0
    error('Initial bubble radius is negative, stop simulation')
end









