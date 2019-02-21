%conserve volume via an eigenfunction

function R = ConserveEigenLinearly(theta,alpha,eigen,r,step)

    %integration
    INT = ([diff(theta),0]+[0,diff(theta)])/2;
% 
%     R = 2/3*pi*INT*((r+alpha*eigen).^3.*sin(theta')) - V0;

    bbb = 2*pi*INT.*(r'.^2.*sin(theta));
    R = bbb*step;
    
    %R = axis_int_gauss_vect((r+alpha*eigen)'.*cos(theta),(r+alpha*eigen)'.*sin(theta)) - V0;
    
end