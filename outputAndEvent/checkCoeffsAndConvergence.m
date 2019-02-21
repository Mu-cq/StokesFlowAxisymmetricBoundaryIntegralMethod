%save results

function status = checkCoeffsAndConvergence(t,xyMode,flag,PARAM,V0)

    %at the first iteration it gives also last point
    if numel(t)==2
        t = t(1);
    end
    
    %save results
    if PARAM.saveDataIte==1
        here = pwd;
        cd(PARAM.res)
        %display('Save data')
        save(PARAM.filenameIte)
        cd(here)
    end
    
    
    status = 0;
    if isempty(t)==1
        status = 1;
    end
    
    %if t is empty (cleaning when exiting) don't compute anything
    if status==0

    %x and y modes
    xMode = xyMode(1:2:end-1);
    yMode = xyMode(2:2:end);
    
    %from modes to grid
    [x,y] = fromModesToGrid(xMode,yMode,PARAM);

    %get smallest coefficient
    xLast = max(abs(xMode(PARAM.dealiasing-1:PARAM.dealiasing)));
    yLast = max(abs(yMode(PARAM.dealiasing-1:PARAM.dealiasing)));
    coeff = max(xLast,yLast);
    
    %condition for stippoing the simulation
    if coeff>PARAM.tresholdCoeffs
        
        disp('Coefficients are too high')
        
        status = 1;
        
    end
    
    %check if dripping occurs
    if PARAM.legendre==1||PARAM.legendre==2
        
        yDrip = y(2:end-1);
        
    elseif PARAM.legendre==0
        
        yDrip = y(2:end-1);
        
    end
    
    if sum(yDrip<0)
        
        disp('Interface escaped the domain')
        
        status = 1;
        
    end
    
    if max(x)>10 || max(y)>4
        
       disp('Droplet is elongating indefitely')
       
       status = 1;
        
    end
    
    if status==0

        %compute residuals
        res = 1;
        if PARAM.checkRes==1

            %check normal velocity
            if PARAM.legendre==1||PARAM.legendre==2
                %compute residulas
                [~,res] = dropLegendreCurvilinearModes(0,xyMode,PARAM);
            elseif PARAM.legendre==0
                [~,res] = dropExtensChebfunCurvilinearModes(0,xyMode,PARAM);
            end

            %check convergence on normal velocity
            if res<PARAM.converge

                disp('Convergence on the residuals has been reached')
                status = 1;

            end
            
            %check shape convergence
            L = max(x)-min(x);
            B = 2*max(y);
            D = (L-B)/(L+B);
            if PARAM.BC==1 && PARAM.convergeShape>0
    
                if PARAM.visc==1
                    load('./steadyState/CaExt')
                    load('./steadyState/DExt')
                elseif PARAM.visc==0
                    load('./steadyState/CaExt0')
                    load('./steadyState/DExt0')
                elseif PARAM.visc==0.02
                    load('./steadyState/CaExt002')
                    load('./steadyState/DExt002')
                elseif PARAM.visc==0.05
                    load('./steadyState/CaExt005')
                    load('./steadyState/DExt005')
                elseif PARAM.visc==0.1
                    load('./steadyState/CaExt01')
                    load('./steadyState/DExt01')
                elseif PARAM.visc==0.5
                    load('./steadyState/CaExt05')
                    load('./steadyState/DExt05')
                elseif PARAM.visc==5
                    load('./steadyState/CaExt5')
                    load('./steadyState/DExt5')
                elseif PARAM.visc==10
                    load('./steadyState/CaExt10')
                    load('./steadyState/DExt10')
                elseif PARAM.visc==0.01
                    load('./steadyState/CaExt001')
                    load('./steadyState/DExt001')
                end
                
                [~,ind] = min(abs(PARAM.Ca-manyCa));
                Dca = manyD(ind);
                if abs(D-Dca)/Dca < PARAM.convergeShape
                    disp('Convergence on the shape has been reached')
                    status  = 1;
                end

            elseif PARAM.BC==2 && PARAM.convergeShape>0

                if abs(D) < PARAM.convergeShape
                    disp('Convergence on the shape has been reached')
                    status  = 1;
                end

            end

        end
        
    %output to screen
    if PARAM.algorithm~=6
        fprintf('T=%f res=%1.16f LastCoeff=%1.16f\n',t,res,coeff);
    end
        
    end
    
    end

end