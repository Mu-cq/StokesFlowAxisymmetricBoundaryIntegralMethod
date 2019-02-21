%find a constant which resize a given drop interface in polar coordinates

function c0 = findC0(theta,r,V0)

    %nonlinear function
    fVolume = @(CCC) ModifyVolumeC0(theta,r,CCC,V0);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    c0 = fsolve(fVolume,1,options);

end