%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [P,Px,Py,P_an] = laplaceST(P,Px,Py,dL)

    %display('ciao')
    for i = 1:numel(dL)
        
        %singularity treatment
        %SINGLE LAYER
        P(1+6*(i-1):6*i,i) = zeros(6,1);
        
        %DOUBLE LAYER
        Px(1+6*(i-1):6*i,i) = zeros(6,1);
        Py(1+6*(i-1):6*i,i) = zeros(6,1);
        
    end
    
    P_an = -1/2/pi*diag(dL-dL.*log(dL/2));

end

