%conserve volume via an eigenfunction

function R = ConserveEigenHalf(theta,alpha,eigen,r,V0)

    %integration
%     INT = ([diff(theta),0]+[0,diff(theta)])/2;
% 
%     R = 2/3*pi*INT*((r+alpha*eigen).^3.*sin(theta')) - V0;

    x = (r+alpha*eigen)'.*cos(theta);
    y = (r+alpha*eigen)'.*sin(theta);
    
    %flip
    x = [x(1:end-1) -flip(x)];
    y = [y(1:end-1) flip(y)];
    
    R = axis_int_gauss_vect(x,y) - V0;
    
end