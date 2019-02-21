%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [P,Px,Py,P_an] = laplaceST(P,Px,Py,dL,x,y,x0,y0)

    %display('ciao')
    for i = 1:numel(dL)
        
        %distance form the singularity
        r = sqrt((x(1+6*(i-1):6*i,i)-x0(i)).^2+(y(1+6*(i-1):6*i,i)-y0(i)).^2);
        
        %singularity treatment
        %SINGLE LAYER
        %P(1+6*(i-1):6*i,i) = P(1+6*(i-1):6*i,i) + log(r)/2/pi/y0(i) + (y(1+6*(i-1):6*i,i)-y0(i)).*log(r)/2/pi/y0(i);
        P(1+6*(i-1):6*i,i) = P(1+6*(i-1):6*i,i) + log(r)/2/pi;
        
        %DOUBLE LAYER
        %Px(1+6*(i-1):6*i,i) = zeros(6,1);
        %Py(1+6*(i-1):6*i,i) = zeros(6,1);
        
    end
    
    %analytycal integration
    P_an = -1/2/pi*diag(dL.*log(dL/2)-dL);

end

