%compute center of mass of a solid of rotation

function xcm = center_mass(x,y)

        %using trapezi rules for integration
        dx = (x(1:end-1)-x(2:end));
                
        num = sum((x(1:end-1).*y(1:end-1).*y(1:end-1)+x(2:end).*y(2:end).*y(2:end)).*dx/2);
        den = sum((y(1:end-1).*y(1:end-1)+y(2:end).*y(2:end)).*dx/2);
        
        xcm = num/den;

return