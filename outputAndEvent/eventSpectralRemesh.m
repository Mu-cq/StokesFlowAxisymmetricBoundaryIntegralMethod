%activate spectral remesh

function [value,isterminal,direction] = eventSpectralRemesh(t,y,PARAM)

%value = round(checkMetricCurvilinear(t,y,PARAM)-1);
value = checkMetricCurvilinear(t,y,PARAM)-1;
isterminal = 1;
direction = 0;