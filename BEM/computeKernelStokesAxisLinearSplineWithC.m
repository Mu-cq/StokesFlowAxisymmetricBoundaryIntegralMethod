%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeKernelStokesAxisLinearSplineWithC(x0,y0,ax,ay,bx,by,cx,cy,dx,dy,visc)

    %number of singularities and elemnts on which the integration is
    %perfromed
    N = numel(x0);
    elem = numel(ax);
    
    %compute
    if visc==1
        
        [GXX,GXY,GYX,GYY] = Stokes2DAxisSPlinesLinearOnlySingle(ax,bx,cx,dx,ay,by,cy,dy,x0,y0);
        
        %reshape because c function gives a vector
        GXX = reshape(GXX,N,elem+1);
        GXY = reshape(GXY,N,elem+1);
        GYX = reshape(GYX,N,elem+1);
        GYY = reshape(GYY,N,elem+1);
        A11 = zeros(N,elem+1);
        A12 = zeros(N,elem+1);
        A21 = zeros(N,elem+1);
        A22 = zeros(N,elem+1);
    
    else
        
        error('It seems that there is a bug in the Double Layer A22, use matlab instead of C to build them matrix')
        [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = Stokes2DAxisSPlinesLinearNoDouble(ax,bx,cx,dx,ay,by,cy,dy,x0,y0);
        
        %reshape because c function gives a vector
        GXX = reshape(GXX,N,elem+1);
        GXY = reshape(GXY,N,elem+1);
        GYX = reshape(GYX,N,elem+1);
        GYY = reshape(GYY,N,elem+1);
        A11 = reshape(A11,N,elem+1);
        A12 = reshape(A12,N,elem+1);
        A21 = reshape(A21,N,elem+1);
        A22 = reshape(A22,N,elem+1);
        
    end
    
    %set NAN to zero
    GXX(isnan(GXX)) = 0;
    GXY(isnan(GXY)) = 0;
    GYX(isnan(GYX)) = 0;
    GYY(isnan(GYY)) = 0;
    A11(isnan(A11)) = 0;
    A12(isnan(A12)) = 0;
    A21(isnan(A21)) = 0;
    A22(isnan(A22)) = 0;

end