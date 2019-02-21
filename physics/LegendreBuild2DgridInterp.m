%compute series of legendre form given data points

function SumP = LegendreBuild2DgridInterp(fff,PPP)
    
    %compute modes corfficiens
    [nS,nT] = size(fff);
    [~,nU] = size(PPP);
    SumP = zeros(size(PPP));
        
      for k = 1:nT
        
            f = fff(:,k);
            f = repmat(f,1,nU);
            
            %legendre function
            P = PPP(k,:);   %legendre polynomia
            P = repmat(P,nS,1);

            %compute serie
            SumP = SumP + f.*P;
        
      end

end