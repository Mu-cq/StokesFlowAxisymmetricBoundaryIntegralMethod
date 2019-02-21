%compute volume flow rate from fixed flux of chemicals

function Q = volumeFlowRate(x,y,Ca,beta,PARAM)

if PARAM.massFlux==0    %impose volume flux
    
   Q = Ca; 
    
elseif PARAM.massFlux==1        %impose mass flux
    
   C0 = 4/3/(3/4/pi)^(1/3);
   V = axis_int_gauss_vect(x,y);
   Q = Ca/(beta+C0*V^(-1/3));
    
end