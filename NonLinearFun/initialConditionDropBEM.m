%choose non linear function, 1 drop, Spectral

function [initialXY,PARAM] = initialConditionDropBEM(PARAM)

%build initial shape
if PARAM.uploadRes==0
    
    disp('Shape is analytical function')
    
    [xInitial,yInitial,PARAM] = buildGeometryPanelsGeneral(PARAM);
    if PARAM.algorithm==2
        error('You have to start from base state obtained from DNS')
    end
    
elseif PARAM.uploadRes==1
    
    disp('Shape is uploaded from previous result')
    
    [xInitial{1},yInitial{1},PARAM] = uploadOneDropOnly(PARAM.ODE,PARAM.TendUP,PARAM.dtUP,PARAM);
    xcm = center_mass(xInitial{1},yInitial{1});
    xInitial{1} = xInitial{1}-xcm;
    
else
    error('Not implemented')
end

%initial condition
initialXY = zeros(2*numel(xInitial{1}),1);
initialXY(1:2:end-1) = xInitial{1};
initialXY(2:2:end) = yInitial{1};

