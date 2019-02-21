%save results

function status = checkCoeffsAndConvergence(t,xyMode,flag,PARAM)

    %x and y modes
    xMode = xyMode(1:2:end-1);
    yMode = xyMode(2:2:end);

    %get smallest coefficient
    if PARAM.volume==1
        xLast = xMode(PARAM.dealiasing-1);
        yLast = yMode(PARAM.dealiasing-1);
        coeff = max(xLast,yLast);
    else
        xLast = xMode(PARAM.dealiasing-1);
        yLast = yMode(PARAM.dealiasing-1);
        coeff = max(xLast,yLast);
    end

    display('Save Data')
    
    status = 0;

end