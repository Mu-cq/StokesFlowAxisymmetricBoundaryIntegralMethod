%green function of Laplace axisymmetric (potential flow source)

function [G,Gx,Gy,NXX,NXY,NYX,NYY] = gf_ax_fs_laplace(x,y,x0,y0,wall1,wall2,visc,p_out,mass_flow)

    K2 = 4*y.*y0./((x-x0).^2+(y-y0).^2);
    [F, E] = ell_int_vect(K);

end
