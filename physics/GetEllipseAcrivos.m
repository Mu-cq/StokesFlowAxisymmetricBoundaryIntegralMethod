% get ellipse shapes from linear theory Acrivos 1973

function D = GetEllipseAcrivos(Ca, lambda)

%CAREFUL BEACUSE MY CAPILLARY NUMBER IS SCALED BY 2
Ca = Ca/2;

%equation parameters
a0 = 5/3/(2*lambda+3);
a1 = -40*(lambda+1)/(2*lambda+3)/(19*lambda+16);
a2 = 10*(4*lambda-9)/7/(2*lambda+3)^2;
a3 = 288*(137*lambda^3+624*lambda^2+741*lambda+248)/(7*(2*lambda+3)^2*(19*lambda+16)^2);

b0 = -360*(lambda+1)/(17*lambda+16)/(10*lambda+11);
b1 = 1/7/(2*lambda+3);
b2 = (16*(-14*lambda^3+207*lambda^2+431*lambda+192))/(21*(2*lambda+3)*(19*lambda+16)*(17*lambda+16)*(10*lambda+11));

F = ((a1+a2*Ca)+sqrt((a1+a2*Ca).^2-4*a0*a3*Ca))./(2*a3*Ca);
D = Ca.*(-4.5*F-399/2/b0*Ca.*F.*(-b1+b2*F))./(1-1.5*Ca.*(F+24/5*Ca.*F.^2)-891/2/b0*Ca.^2.*F.*(-b1+b2*F));

end