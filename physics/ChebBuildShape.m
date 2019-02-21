%compute series of legendre form given data points

function SumP = ChebBuildShape(theta,f,symmetric)
    
    %compute modes corfficiens
    %SumP = f(1)/2;
    modes = numel(f);
    
    SumP = -0.5*f(1);
    
    for i = 1:modes

            %I take only symmetric modes if the shape is symmetric
            temp = i-1;
            if symmetric==1
                polN = 2*temp;
            elseif symmetric==0
                polN = temp;
            end

            %chebishev
            T = cos(polN*theta);
            
            %compute serie
            SumP = SumP + f(polN+1)*T;
        
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