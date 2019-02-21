%replace interface point too close to walls

function df = repulsiveStress(PARAM,distDrop)

    %display('displace')
    active = distDrop<PARAM.repulsiveOn;
    
    %compute repulsive force
    df = PARAM.coeffRepulsive*repulsiveForce(distDrop,PARAM).*active;

end