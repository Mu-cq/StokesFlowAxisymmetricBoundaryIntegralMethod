% get ellipse shapes from linear theory Acrivos 1973

function R = GetParabolaAcrivos(alpha,t,Ca,V0)

x = (-t+1/2)*alpha;
a = max(x);
y = ParabolaAcrivos(a,x,Ca);

%compute volume
V = VolumeCurvilinearAxisSpectral(x,y,PARAM);

R = abs(V-V0)/V0;