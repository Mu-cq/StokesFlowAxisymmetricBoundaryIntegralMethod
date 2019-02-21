%replace interface point too close to walls

function df = disjoiningPressure(PARAM,distDrop)

    %compute repulsive force
    df = PARAM.hamacker./distDrop.^3;

end