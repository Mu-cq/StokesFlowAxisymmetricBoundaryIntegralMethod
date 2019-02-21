%Choose bubble position and size

function [PARAM,value] = chooseBubblePositionAndSize(PARAM,beta,Hcc,nPerLenght,bubbleNucleates)

%PARAM here
value = 0;
PARAMhere = PARAM;

%modify parameters for cone without bubble
PARAMhere.panels = 4;
PARAMhere.n = PARAM.n(1:4);

%build geometry
[x,y,PARAMhere] = buildGeometryPanelsGeneral(PARAMhere);

%solve Laplace equation
yLaplace = BEM_Laplace(x,y,PARAMhere);

%compute concentration along the axis
X = linspace(-50,60,500);
Y = zeros(1,500);
[~,~,PHIfield] = computeConcentrationField(X,Y,x,y,yLaplace,PARAMhere,0,1);

%find maximum
[maxPHI,ind] = max(PHIfield);
xcmBubble = X(ind);
PARAM.maxConcAxis = bubbleNucleates;

if maxPHI<bubbleNucleates
    
   display('Concentration is to low for bubble nucleation to happpen')
   value = 1;
    
end

%compute minimum radius from Henry's law
%rMin = 2/Hcc/(bubbleNucleates-beta/Hcc);
rMin = 2*Hcc/(maxPHI-beta*Hcc);
rMin = rMin+0.1;

%check if bubble radius is too large
if rMin>PARAM.maxRadiusNucleation
    display(['Initial radius it too large, r=' num2str(rMin)])
    value = 1;
end

if rMin<PARAM.smallestRadius
    rMin = PARAM.smallestRadius;
    display(['Minimum bubble size allowed is r=' num2str(PARAM.smallestRadius)])
end

%check if bubble radius is too large
if rMin==0
    display('Initial bubble radius is negative, stop simulation')
    value = 1;
end

%solve Laplace equation with bubble
PARAMbubble = PARAM;
PARAMbubble.rArc(5) = rMin; PARAMbubble.x0_Circle(5) = xcmBubble;
PARAMbubble.n(5) = round(rMin*pi*20);
if PARAMbubble.n(5)<10
    PARAMbubble.n(5) = 10;
elseif PARAMbubble.n(5)>100
    PARAMbubble.n(5) = 100;
end
[xBubble,yBubble,PARAMbubble] = buildGeometryPanelsGeneral(PARAMbubble);
PARAMbubble.concBC{5} = Hcc*(beta+2/rMin);
yLaplaceBubble = BEM_Laplace(xBubble,yBubble,PARAMbubble);

%compute mass flow rate to the bubble
Qmass = computeFlowRate(xBubble{5},yBubble{5},yLaplaceBubble(sum(PARAMbubble.n(1:4))+1:sum(PARAMbubble.n)),5,PARAMbubble);

if Qmass<0
    
   display('Bubble is shrinking at the first iteration')
   value = 1;
    
end

if value==0

    display(['Bubble nucleates at z=' num2str(xcmBubble) ' with radius r=' num2str(rMin)])

end

%push to PARAM
PARAM.x0_Circle(5) = xcmBubble;
PARAM.rArc(5) = rMin;
PARAM.n(5) = round(rMin*pi*nPerLenght);
if PARAM.n(5)<10
    PARAM.n(5) = 10;
end
PARAM.maxElem(5) = 1/nPerLenght;




