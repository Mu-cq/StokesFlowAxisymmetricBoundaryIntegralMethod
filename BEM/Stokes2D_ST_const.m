%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [SXX,SYY,QXXX,QXXY,QXYY,QYYY,GXX_an,GYY_an] =...
    Stokes2D_ST_const(SXX,SYY,QXXX,QXXY,QXYY,QYYY,walls,globalX,globalY,X0,Y0,dL,elements)

    for i = 1:walls
        
        %singularity treatment
        %SINGLE LAYER
        SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i)+log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))';
        SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i)+log(sqrt((globalX(1+6*(i-1):6*i)-X0(i)).^2+(globalY(1+6*(i-1):6*i)-Y0(i)).^2))';
        
        %DOUBLE LAYER
        QXXX(1+6*(i-1):6*i,i) = zeros(6,1);
        QXXY(1+6*(i-1):6*i,i) = zeros(6,1);
        QXYY(1+6*(i-1):6*i,i) = zeros(6,1);
        QYYY(1+6*(i-1):6*i,i) = zeros(6,1);
        
    end
    
    GXX_an = [diag(-dL(1:walls).*log(dL(1:walls)/2)+dL(1:walls)); zeros(elements-walls,walls)];
    GYY_an = GXX_an;

end

