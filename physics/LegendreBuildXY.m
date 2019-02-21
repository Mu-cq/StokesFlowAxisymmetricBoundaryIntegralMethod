%compute series of legendre form given data points

function SumP = LegendreBuildXY(f,PPP)
    
    %compute modes corfficiens
    SumP = 0;
    modes = numel(f);
    for i = 1:modes

            %legendre function
            P = PPP(i,:)';   %legendre polynomia

            %compute serie
            SumP = SumP + f(i)*P;
        
    end

end