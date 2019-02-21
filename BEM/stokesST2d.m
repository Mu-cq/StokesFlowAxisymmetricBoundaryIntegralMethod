%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [SXX,SYY,QXXX,QXXY,QXYY,QYYY,SXXan,SYYan] = stokesST2d(SXX,SYY,QXXX,QXXY,QXYY,QYYY,dL,x,y,x0,y0)

    %display('ciao')
    for i = 1:numel(dL)
        
        %distance form the singularity
        dx = x(1+6*(i-1):6*i,i)-x0(i);
        dy = y(1+6*(i-1):6*i,i)-y0(i);
        r = sqrt((x(1+6*(i-1):6*i,i)-x0(i)).^2+(y(1+6*(i-1):6*i,i)-y0(i)).^2);
        
        %singularity treatment
        %SINGLE LAYER
        SXX(1+6*(i-1):6*i,i) = dx.^2./r.^2;
        SYY(1+6*(i-1):6*i,i) = dy.^2./r.^2;
        
        %DOUBLE LAYER
        QXXX(1+6*(i-1):6*i,i) = zeros(6,1);
        QXXY(1+6*(i-1):6*i,i) = zeros(6,1);
        QXYY(1+6*(i-1):6*i,i) = zeros(6,1);
        QYYY(1+6*(i-1):6*i,i) = zeros(6,1);
        
    end
    
    %analytycal integration single layer
    SXXan = -diag(dL.*log(dL/2)-dL);
    SYYan = -diag(dL.*log(dL/2)-dL);

end

