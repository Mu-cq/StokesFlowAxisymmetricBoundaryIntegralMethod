%given a certain amplitude of area variation norm perturbation, find the
%corresponding value of prolate/oblate perturabtion with the constrain of
%the volume being 4*pi/3

function ell = elliptic_perturb(a_perturb,pro_obl)

    if pro_obl==1
        %formula for prolates (APPROXIMATED FORMULA FOR THE AREA!!!)
        %syms x
        %a = fsolve(@(x) 2*pi/x*(1+(2*x^3-x^6)/(1-x^3)*atanh(1-x^3))-8*pi-a_perturb,0.5);
        a = fsolve(@(x) 4*pi*((2*x^0.8+x^(-1.6))/3)^(1/1.6)-4*pi-a_perturb,1.1);
        b = sqrt(1/a);
        ell = abs((a-b)/(a+b));
    elseif pro_obl==2
        %formulas for oblates
        
        
    end

end