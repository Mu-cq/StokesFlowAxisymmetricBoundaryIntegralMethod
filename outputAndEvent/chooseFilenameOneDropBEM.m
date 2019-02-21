%choose filename

function filename = chooseFilenameOneDropBEM(PARAM,maxDT,Tend,A0,ODE,nameShort)

%PARAM.filename = ['DropSpectral_ODE=' num2str(PARAM.ODE) '_Legendre=' num2str(PARAM.legendre) '_BC=' num2str(PARAM.BC) '_Ca=' num2str(PARAM.Ca) '_visc=' num2str(PARAM.visc) '_n=' num2str(PARAM.dealiasing) '_D=' num2str(PARAM.D) '_maxDT=' num2str(maxDT) '_VolCorr=' num2str(PARAM.volume) '_Tend=' num2str(Tend) '_CPUs=' num2str(cpu) '.mat'];
if isempty(PARAM.Ca)
    
    if PARAM.algorithm==1
        filename = [nameShort 'oneDropBEM_ODE=' num2str(ODE) '_n=' num2str(PARAM.n) '_BC=' num2str(PARAM.typeBCstokes) '_Bo=' num2str(PARAM.Bond) '_visc=' num2str(PARAM.visc) '_D=' num2str(PARAM.D) '_maxDT=' num2str(maxDT) '_Tend=' num2str(Tend) '.mat'];
    elseif PARAM.algorithm==2
        filename = [nameShort 'minimalSeedBEM_n=' num2str(PARAM.n) '_Bo=' num2str(PARAM.Bond) '_visc=' num2str(PARAM.visc) '_BC=' num2str(PARAM.typeBCstokes) '_T=' num2str(Tend) '_A0=' num2str(A0) '.mat'];
    elseif PARAM.algorithm==3
        filename = ['edgeTrackingDropBEM_edgeLoop=' num2str(PARAM.edgeLoop) '_deltaEdge=' num2str(PARAM.deltaEdge) '_ODE=' num2str(PARAM.ODE) '_BC=' num2str(PARAM.typeBCstokes) '_Bo=' num2str(PARAM.Bond) '_visc=' num2str(PARAM.visc) '_n=' num2str(PARAM.n) '_maxDT=' num2str(maxDT) '.mat'];
    else
        error('Not implemented')
    end
    
elseif isempty(PARAM.Bond)
    
    if PARAM.algorithm==1
        filename = [nameShort 'oneDropBEM_ODE=' num2str(ODE) '_n=' num2str(PARAM.n) '_BC=' num2str(PARAM.typeBCstokes) '_Ca=' num2str(PARAM.Ca) '_visc=' num2str(PARAM.visc) '_D=' num2str(PARAM.D) '_maxDT=' num2str(maxDT) '_Tend=' num2str(Tend) '.mat'];
    elseif PARAM.algorithm==2
        filename = [nameShort 'minimalSeedBEM_n=' num2str(PARAM.n) '_Ca=' num2str(PARAM.Ca) '_visc=' num2str(PARAM.visc) '_BC=' num2str(PARAM.typeBCstokes) '_T=' num2str(Tend) '_A0=' num2str(A0) '.mat'];
    else
        error('Not implemented')
    end
    
    
else
    
    error('The name for the output has not been chosen')
    
end