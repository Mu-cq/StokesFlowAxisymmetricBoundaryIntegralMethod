%compute series of legendre form given data points

function SumP = LegendreBuild2DgridCosTheta(fff,PPP)
    
    %compute modes corfficiens
    [nS,nT] = size(fff);
    %[~,nV] = size(PPP);
    SumP = zeros(size(fff));
        
      for k = 1:nT
        
            f = fff(:,k);
            f = repmat(f,1,nT);
            
            %legendre function
            P = PPP(k,:);   %legendre polynomia
            P = repmat(P,nS,1);

            %compute serie
            SumP = SumP + f.*P;
        
      end

end