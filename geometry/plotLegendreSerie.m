%plot shape from legendre modes

function [xNew,yNew] = plotLegendreSerie(x,y,fMode,symmetric)

xCM = center_mass(x,y);
x = x-xCM;

% from cartesian to polar coordinates
theta = atan(y./x);
theta = theta + pi*(theta<0);

modes = numel(fMode);

%compute new shape form legendre serie
r = 0;
for i = 1:modes
        
        %I take only symmetric modes if the shape is symmetric
        temp = i-1;
        if symmetric==1
            polN = 2*temp;
        elseif symmetric==0
            polN = temp;
        end

        %legendre function
        PPP = legendre(polN,cos(theta));
        P = PPP(1,:)';   %legendre polynomia
        
        r = r + fMode(i)*P;
        
end

xNew = r.*cos(theta) + xCM;
yNew = r.*sin(theta);