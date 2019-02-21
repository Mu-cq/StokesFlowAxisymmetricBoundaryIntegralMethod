%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolume(a,b,rho,V0)

    n = numel(a);
    x = [a(1:n/2) 0 a(n/2+1:end)];
    y = [b(1:n/2) rho b(n/2+1:end)];
    
    R = axis_int_gauss_vect(x,y)-V0;

end