%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [SXX,SYY,QXXX,QXXY,QXYY,QYYY,GXX_an,GYY_an] =...
                Stokes2D_ST_linear(SXX,SYY,QXXX,QXXY,QXYY,QYYY,walls,globalX,globalY,X0,Y0,dL,el)
    
    %display('ciao lin sing')

    % singularity treatment lin elem
    for i = walls+1:el-1
        
       dx1 = globalX(1+6*(i-1):6*i)-X0(i);
       dy1 = globalY(1+6*(i-1):6*i)-Y0(i);
       dx2 = globalX(1+6*(i-1):6*i)-X0(i+1);
       dy2 = globalY(1+6*(i-1):6*i)-Y0(i+1);
       
       %SINGLE LAYER
       r2 = dx1.^2+dy1.^2;
       SXX(1+6*(i-1):6*i,i) = (dx1.*dx1)./r2;
       SYY(1+6*(i-1):6*i,i) = (dy1.*dy1)./r2;
       
       r2 = dx2.^2+dy2.^2;
       SXX(1+6*(i-1):6*i,i+1) = (dx2.*dx2)./r2;
       SYY(1+6*(i-1):6*i,i+1) = (dy2.*dy2)./r2;
       
       QXXX(1+6*(i-1):6*i,i) = zeros(6,1);
       QXXY(1+6*(i-1):6*i,i) = zeros(6,1);
       QXYY(1+6*(i-1):6*i,i) = zeros(6,1);
       QYYY(1+6*(i-1):6*i,i) = zeros(6,1);
       QXXX(1+6*(i-1):6*i,i+1) = zeros(6,1);
       QXXY(1+6*(i-1):6*i,i+1) = zeros(6,1);
       QXYY(1+6*(i-1):6*i,i+1) = zeros(6,1);
       QYYY(1+6*(i-1):6*i,i+1) = zeros(6,1);
        
    end
    
    %because last element touches first node PRIODICITY CONDITION
    %SINGLE LAYER
    dxel = globalX(1+6*(el-1):6*el)-X0(el);
    dyel = globalY(1+6*(el-1):6*el)-Y0(el);
    
    r2 = dxel.^2+dyel.^2;
    SXX(1+6*(el-1):6*el,el) = (dxel.*dxel)./r2;
    SYY(1+6*(el-1):6*el,el) = (dyel.*dyel)./r2;
    QXXX(1+6*(el-1):6*el,el) = zeros(6,1);
    QXXY(1+6*(el-1):6*el,el) = zeros(6,1);
    QXYY(1+6*(el-1):6*el,el) = zeros(6,1);
    QYYY(1+6*(el-1):6*el,el) = zeros(6,1);
    
    dxel1 = globalX(1+6*(el-1):6*el)-X0(walls+1);
    dyel1 = globalY(1+6*(el-1):6*el)-Y0(walls+1);
    
    r2 = dxel1.^2+dyel1.^2;
    SXX(1+6*(el-1):6*el,walls+1) = (dxel1.*dxel1)./r2;
    SYY(1+6*(el-1):6*el,walls+1) = (dyel1.*dyel1)./r2;
    QXXX(1+6*(el-1):6*el,walls+1) = zeros(6,1);
    QXXY(1+6*(el-1):6*el,walls+1) = zeros(6,1);
    QXYY(1+6*(el-1):6*el,walls+1) = zeros(6,1);
    QYYY(1+6*(el-1):6*el,walls+1) = zeros(6,1);

    %adding analytical term following the notation in the notes
    A = 3/4*dL-0.5*dL.*(log(dL));
    B = 1/4*dL-0.5*dL.*(log(dL));
    
    %build matrix
    GXX_an = diag(A + [0 A(1:end-1)]) + diag(B(1:end-1),1) + diag(B(1:end-1),-1);
    %periodicity
    GXX_an(1,1) = GXX_an(1,1) + A(end);
    GXX_an(1,end) = B(end); GXX_an(end,1) = B(end);
    
    GYY_an = GXX_an;
    

end

