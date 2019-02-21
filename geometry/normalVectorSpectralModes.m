%compute normal vector speectrally

function [nx,ny] = normalVectorSpectralModes(xMode,yMode,PARAM)

    %derivatives
    D1 = PARAM.D1;
    if PARAM.legendre==1
        PPP = PARAM.PPP;
    end

    %transform in grid points
    if PARAM.legendre==0
        x = chebvals2chebcoeffs(xMode);
        y = chebvals2chebcoeffs(yMode);
    elseif PARAM.legendre==1
        x = LegendreBuildXY(xMode,PPP);
        y = LegendreBuildXY(yMode,PPP);
    end

    %compute geomtrical derivaties
    xp = D1*x;    yp = D1*y;

    %compute normal vector
    h = (xp.^2+yp.^2).^(0.5);
    nx = yp./h;
    ny = -xp./h;

    %compute modes
    if PARAM.legendre==0
        nx = chebfun(nx,[0 1]);
        ny = chebfun(ny,[0 1]);
        nx = chebcoeffs(nx);
        ny = chebcoeffs(ny);
    elseif PARAM.legendre==1
        nx = LegendreSerieSpectralXY(nx,PPP,PARAM);
        ny = LegendreSerieSpectralXY(ny,PPP,PARAM);
    end
    
end