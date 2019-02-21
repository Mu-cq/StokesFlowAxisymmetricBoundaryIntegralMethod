%replace interface point too close to walls

function df = disjoiningPressureDoubleLayer(PARAM,distDrop)

    delta = PARAM.repulsiveOn;
    r = distDrop;

    %compute only when surfaces are close
    activate = distDrop<=delta;

    %compute repulsive force
    df = PARAM.coeffRepulsive*(exp(delta-r)-1)./(exp(delta)-1).*activate;

end