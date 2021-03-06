%compute series of legendre form given data points

function f = LegendreSerie(x,y,modes,symmetric,PARAM)

    %weight = PARAM.weightSpline;

    %plaxe x in the center of mass
    xCM = center_mass(x,y);
    x = x-xCM;

    % from cartesian to polar coordinates
    r = sqrt(x.^2+y.^2);
    theta = atan(y./x);
    theta = theta + pi*(theta<0);

    %initialize variables
    f = zeros(1,numel(modes));
    
    %compute modes corfficiens
    SumP = 0;
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

            %compute coeff
            int = r.*P.*sin(theta);
            f(i) = (2*polN+1)/2*trapz(theta,int);
            %f(i) = (2*polN+1)/2*weight*int;

            %compute serie
            SumP = SumP + f(i)*P;
        
    end
    
%     figure
%     plot(theta,r)
%     hold on
%     plot(theta,SumP,'--')
%     grid on
%     hold off
    
    %Vserie = 2/3*INT*(r.^3.*sin(theta));
    %err = abs(Vserie-V0)/V0;

end