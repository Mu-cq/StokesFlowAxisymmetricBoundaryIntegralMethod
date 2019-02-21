%save results

function status = saveResults(t,y,flag,PARAM)

    display('Save Data')
    
    if numel(t)==2
        t = t(1);
    end

    %save file
    filename = ['t=' num2str(t) '.mat'];
    cd(PARAM.res)
    save(filename,'y')
    cd(PARAM.this_folder)
    status = 0;

end