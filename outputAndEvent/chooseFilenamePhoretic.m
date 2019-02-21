%choose filename

function filename = chooseFilenamePhoretic(PARAM)

%PARAM.filename = ['DropSpectral_ODE=' num2str(PARAM.ODE) '_Legendre=' num2str(PARAM.legendre) '_BC=' num2str(PARAM.BC) '_Ca=' num2str(PARAM.Ca) '_visc=' num2str(PARAM.visc) '_n=' num2str(PARAM.dealiasing) '_D=' num2str(PARAM.D) '_maxDT=' num2str(maxDT) '_VolCorr=' num2str(PARAM.volume) '_Tend=' num2str(Tend) '_CPUs=' num2str(cpu) '.mat'];
if PARAM.algorithm==1
    filename = ['DropSpectral_ODE=' num2str(PARAM.ODE) '_Legendre=' num2str(PARAM.legendre) '_BC=' num2str(PARAM.BC) '_Ca=' num2str(PARAM.Ca) '_visc=' num2str(PARAM.visc) '_n=' num2str(PARAM.dealiasing) '_D=' num2str(PARAM.D) '_maxDT=' num2str(maxDT) '_VolCorr=' num2str(PARAM.volume) '_Tend=' num2str(Tend) '.mat'];
elseif PARAM.algorithm==2
    filename = ['edgeTrackingDrop_edgeLoop=' num2str(PARAM.edgeLoop) '_deltaEdge=' num2str(PARAM.deltaEdge) '_ODE=' num2str(PARAM.ODE) '_Legendre=' num2str(PARAM.legendre) '_BC=' num2str(PARAM.BC) '_Ca=' num2str(PARAM.Ca) '_visc=' num2str(PARAM.visc) '_n=' num2str(PARAM.dealiasing) '_maxDT=' num2str(maxDT) '_VolCorr=' num2str(PARAM.volume) '.mat'];
elseif PARAM.algorithm==3
    filename = ['newtonMethodSpectralXYmodes_n=' num2str(PARAM.dealiasing) '_Ca=' num2str(PARAM.Ca) '_visc=' num2str(PARAM.visc) '_Legendre=' num2str(PARAM.legendre) '_BC=' num2str(PARAM.BC) '.mat'];
elseif PARAM.algorithm==4
    filename = ['newtonMethodSpectralXYmodesCont_n=' num2str(PARAM.dealiasing) '_CaUp=' num2str(PARAM.CaBreakUp) '_CaDown='  num2str(PARAM.CaBreakDown) '_visc=' num2str(PARAM.visc) '_Legendre=' num2str(PARAM.legendre) '_BC=' num2str(PARAM.BC) '.mat'];
elseif PARAM.algorithm==5
    filename = [];
else
    error('The name for the output has not been chosen')
end