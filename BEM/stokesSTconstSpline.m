%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [SXX,SYY,anX,anY] = stokesSTconstSpline(SXX,SYY,eta,beta,h,h0,x0,y0,GPt)

    %singularity treatment
    for i = 1:numel(h0)
        
       rangeHere = 1+6*(i-1):6*i;
        
       %distance form the singularity
       r1 = sqrt((eta(rangeHere)-x0(i)).^2+(beta(rangeHere)-y0(i)).^2);
       
       %single layer
       SXX(rangeHere,i) = SXX(rangeHere,i) + (2*log(r1).*h(rangeHere)...
           - 2*(log(r1./GPt).*h(rangeHere) + log(GPt).*(h(rangeHere)-h0(i))))';
       
       SYY(rangeHere,i) = SYY(rangeHere,i) + (2*log(r1).*h(rangeHere)...
           - 2*(log(r1./GPt).*h(rangeHere) + log(GPt).*(h(rangeHere)-h0(i))))';
       
    end
    
    %analytical integration
    anX = diag(2*h0);
    anY = anX;

end

