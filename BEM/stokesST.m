%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [SXX,SYY,SXXan,SYYan] = stokesST(SXX,SYY,dL,x,y,x0,y0)

    %display('ciao')
    for i = 1:numel(dL)
        
        %distance form the singularity
        r = sqrt((x(1+6*(i-1):6*i,i)-x0(i)).^2+(y(1+6*(i-1):6*i,i)-y0(i)).^2);
        
        %singularity treatment
        %SINGLE LAYER
        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + 2*log(r);
        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + 2*log(r);
        
    end
    
    %analytycal integration single layer
    SXXan = -2*diag(dL.*log(dL/2)-dL);
    SYYan = -2*diag(dL.*log(dL/2)-dL);

end

