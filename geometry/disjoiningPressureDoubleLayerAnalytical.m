%replace interface point too close to walls

function [dfX12,dfY12] = disjoiningPressureDoubleLayerAnalytical(a1,b1,a2,b2,PARAM)

    %compute distances and direction
    [dist,nx,ny,dl1,Y1] = geometryRepulsiveDoubleLayerAnalytical(a1,b1,a2,b2,PARAM.theta);

    %critical distance
    delta = PARAM.repulsiveOn;
    r = dist;

    %compute only when surfaces are close
    activate = r<=delta;

    %compute intensity repulsive force
    df = PARAM.coeffRepulsive*(exp(delta-r)-1)./(exp(delta)-1).*activate;
    dfX = -df.*nx;
    dfY = -df.*ny;
    
    %force from 1 to 2
    dfX12 = 2*pi*sum(dfX.*dl1.*Y1,2);
    dfY12 = 2*pi*sum(dfY.*dl1.*Y1,2);

end