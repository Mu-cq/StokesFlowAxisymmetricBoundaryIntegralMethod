%compute series of legendre form given data points

function SumP = FourierBuild1D(an,bn,s)
        
    SumP = an(1)/2;
   for k = 2:numel(an)
        
       n = k-1;
       SumP = SumP + an(k)*cos(n*s) + bn(k)*sin(n*s);
        
   end

end