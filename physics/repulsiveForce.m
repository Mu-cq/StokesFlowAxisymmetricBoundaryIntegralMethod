%compute repulsive forces, exponentail law

function df = repulsiveForce(x,PARAM)

    %activation distance
    repuls = PARAM.repulsiveOn;
    
    %compute repulsive stress
    df = 1./(exp(x/repuls)-1)-1/(exp(1)-1);

end