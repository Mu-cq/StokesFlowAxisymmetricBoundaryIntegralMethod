%compute series of legendre form given data points

function SumP = LegendreBuildShapeChebfun(theta,f,symmetric)
    
    %compute modes corfficiens
    SumP = 0;
    modes = numel(f);
    for i = 1:modes

            %i = modes-k+1;
        
            %I take only symmetric modes if the shape is symmetric
            temp = i-1;
            if symmetric==1
                polN = 2*temp;
            elseif symmetric==0
                polN = temp;
            end

            %legendre function
            PPP = legendre(polN,cos(theta));
            P = PPP(1,:)';   %legendre polynomial

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