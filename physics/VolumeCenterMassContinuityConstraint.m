%compute velocity normal to the interface in an extensional flow

function [DisEquality,Equality] = VolumeCenterMassContinuityConstraint(theta,r,Vfinal,xcmFinal,aWall,bWall,PARAM)
    
    %droplet coordinates
    aDrop = r'.*cos(theta);
    bDrop = r'.*sin(theta);

    %parameters
    n = PARAM.n;   m = PARAM.m;    j = PARAM.j;

    %trapezi integral
    INT = ([diff(theta),0]+[0,diff(theta)])/2;
    
    %volume constraint
    DisEquality = [];
    EqualityVol = 2/3*pi*INT*(r.^3.*sin(theta'))-Vfinal;
    EqualityCM = 3/4 * (INT*(r.^4.*sin(theta').*cos(theta')))/(INT*(r.^3.*sin(theta'))) -xcmFinal;
    
    %compute continuity as integral equation
    a = [aWall aDrop];  b = [bWall bDrop];
    [y,~,~,~,N] = BEM_drop_channel(a,b,PARAM);  %solution
    ux = y(2*n+2*m+2*j+1:2:end-1);  uy = y(2*n+2*m+2*j+2:2:end);    %velocities
    u = N(1,:)'.*ux + N(2,:)'.*uy;
    Vdrop = DropVelocityAxis(aDrop,bDrop,u);
    u = N(1,:)'.*(ux-Vdrop) + N(2,:)'.*uy;
    EqualityContinuity = int_axis_spline_symmetric(aDrop,bDrop,u);
    
    Equality = [EqualityVol EqualityCM EqualityContinuity];

end