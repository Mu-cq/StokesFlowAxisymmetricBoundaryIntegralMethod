%compute series of legendre form given data points

function SumP = FourierBuild2Dgrid(an,bn,s)
    
    %compute modes corfficiens
    [nS,nT] = size(an);
    %[nS,nT] = size(an);
    SumP = zeros(size(an));
    
    n = 0:nS-1;
    nnn = repmat(n',1,nS);
    manyCos = cos(repmat(s',nS,1).*nnn);
    manyCos(1,:) = manyCos(1,:)/2;
    manySin = sin(repmat(s',nS,1).*nnn);
        
      for k = 1:nS
        
            anHere = an(k,:);
            anHere = repmat(anHere,nS,1);
            bnHere = bn(k,:);
            bnHere = repmat(bnHere,nS,1);
            
            %projection function
            cosHere = manyCos(k,:);
            cosHere = repmat(cosHere',1,nT);
            sinHere = manySin(k,:);
            sinHere = repmat(sinHere',1,nT);

            %compute serie
            SumP = SumP + anHere.*cosHere + bnHere.*sinHere;
        
      end

end