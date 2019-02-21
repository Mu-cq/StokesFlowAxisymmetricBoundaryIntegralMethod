%determine the smallest time scale of the problem

function T = time_scale(Ca,R,U,visc)

    T_visc = R/U*(1+visc);   %viscous time scale
    T_surftens = T_visc*Ca; %surface tension time scale
    
    T = min(T_surftens,T_visc); %select the BIGGER one

end