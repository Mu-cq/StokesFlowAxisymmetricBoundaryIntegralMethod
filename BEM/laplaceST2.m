%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [P,PxPy,P_an,PxPy_an] = laplaceST2(P,PxPy,dL,x,y,x0,y0,cosTheta)

    %display('ciao')
    for i = 1:numel(dL)
        
        %distance form the singularity
        r = sqrt((x(1+6*(i-1):6*i,i)-x0(i)).^2+(y(1+6*(i-1):6*i,i)-y0(i)).^2);
        
        %singularity treatment
        %SINGLE LAYER
        P(1+6*(i-1):6*i,i) = P(1+6*(i-1):6*i,i) + log(r)/2/pi;
        
        %DOUBLE LAYER
        PxPy(1+6*(i-1):6*i,i) = PxPy(1+6*(i-1):6*i,i) - cosTheta(i)/(4*pi*y0(i))*log(r);
        
    end
    
    %analytycal integration single layer
    P_an = -1/2/pi*diag(dL.*log(dL/2)-dL);
    
    %analytycal integration double layer
    PxPy_an = 1/4/pi*diag((dL.*log(dL/2)-dL)./y0'.*cosTheta);

end

