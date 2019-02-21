%choose non linear function, 1 drop, Spectral

function f = functionDropSpectral(PARAM)

if PARAM.legendre==1||PARAM.legendre==2
        f = @(t,xyMode) dropLegendreCurvilinearModes(t,xyMode,PARAM);
elseif PARAM.legendre==0
        f = @(t,xyMode) dropExtensChebfunCurvilinearModes(t,xyMode,PARAM);
end