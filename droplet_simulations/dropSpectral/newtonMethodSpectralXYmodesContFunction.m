%compute interface shape for droplet in extensional flow usinf Newton
%method

function newtonMethodSpectralXYmodesContFunction(PARAM)

%current directory
here = pwd;

%print to screen
printToScreen(PARAM);

%local variables
delta = PARAM.delta;                                               % continaution method step
CaBreakUp = PARAM.CaBreakUp;                                           % break when Ca is larger
CaBreakDown = PARAM.CaBreakDown;                                          % break when Ca is smaller
changeSystemOfCoord = PARAM.changeSystemOfCoord;    howOften = PARAM.howOften;                   % change system of coordinates every continuation iteration
switchBranch = PARAM.switchBranch;   deltaSwitch  = PARAM.deltaSwitch;    deltaDir = PARAM.deltaDir;   % for branch switching
isSing = PARAM.isSing;      checkOnlyOnePitchfork = PARAM.checkOnlyOnePitchfork;

%often used
if PARAM.legendre==0
    PPP = PARAM.TTT;
elseif PARAM.legendre==1||PARAM.legendre==2
    PPP = PARAM.PPP;
end
col = 'krgmby';
col = repmat(col,1,10);
PARAM.ODE = 0;

%initial condition
[xyMode,x,y,newCa] = initialConditionDrop(PARAM);
if PARAM.uploadShape==1
    PARAM.Ca = newCa-delta;
else
    PARAM.Ca = PARAM.Ca-delta;
end

%coordinates for plotting
aIN = x;  bIN = y;

%dealiasing
xBase = xyMode(1:2:end-1);
yBase = xyMode(2:2:end);
xBase = xBase(1:PARAM.dealiasing);
yBase = yBase(1:PARAM.dealiasing);

%final volume that I want
V0 = PARAM.V0;
xcm = CenterMassCurvAxisSpectral(x,y,PARAM);

%initialization for saving data
guessLoops = 5e3;
manyCa = zeros(guessLoops,1);
manyD = zeros(guessLoops,1);
manyA = zeros(PARAM.n+1,guessLoops);
manyB = zeros(PARAM.n+1,guessLoops);
minEig = 1;
unstEigen = zeros(guessLoops,1);
manyXbase = zeros(PARAM.dealiasing,guessLoops);
manyYbase = zeros(PARAM.dealiasing,guessLoops);

%initialize for Newton Method
if PARAM.dropFrame==0
    perturb = zeros(PARAM.dealiasing-1,1);
    dir = [zeros(PARAM.dealiasing-1,1); 1]; % initial direction of continuation
    R = zeros(PARAM.dealiasing,1);
    manyPerturb = zeros(PARAM.dealiasing-1,guessLoops);
elseif PARAM.dropFrame==1
    perturb = zeros(PARAM.dealiasing-2,1);
    dir = [zeros(PARAM.dealiasing-2,1); 1]; % initial direction of continuation
    R = zeros(PARAM.dealiasing-1,1);
    manyPerturb = zeros(PARAM.dealiasing-2,guessLoops);
else
    error('Choose if the drop is in the frame of reference of the lab is the co-moving one')
end
sol = [perturb; PARAM.Ca];

%nonlinear equation
fNonlinear = @(unk) NormalVelocitySpectralVolumeXYmodesCont(unk,xBase,yBase,V0,xcm,PARAM,sol(end));
[~,nx,ny,xBaseGrid,yBaseGrid] = fNonlinear(sol(1:end-1));

%check initial shape
if min(ny)<-1e-8
    warning('Normal vector of the initial shape is not accurate')
end

%for saving
%filename = ['newtonMethodSpectralXYmodesCont_n=' num2str(PARAM.dealiasing) '_CaUp=' num2str(CaBreakUp) '_CaDown='  num2str(CaBreakDown) '_visc=' num2str(PARAM.visc) '.mat'];
PARAM.filename = ['newtonMethodSpectralXYmodesCont_n=' num2str(PARAM.dealiasing) '_CaUp=' num2str(CaBreakUp) '_CaDown='  num2str(CaBreakDown) '_visc=' num2str(PARAM.visc) '_Legendre=' num2str(PARAM.legendre) '_BC=' num2str(PARAM.BC) '.mat'];

if PARAM.dropFrame==1
    error('Bug in remeshing when droplet move in its own frame')
end

%continuation
display('Continuation loop')
countCont = 1;
beforeSing = 0;
while (PARAM.Ca<CaBreakUp && PARAM.Ca>CaBreakDown) || countCont==1
    
    if changeSystemOfCoord==1 && sum(countCont==0:howOften:2000)
        
        disp('Update base state shape')
        
        %present solution
        perturb = sol(1:end-1); PARAM.Ca = sol(end);
        if PARAM.dropFrame==0
            x = xBaseGrid + PARAM.CutStep*nx.*sum(repmat(perturb,1,numel(x)).*PPP(2:PARAM.dealiasing,:))';
            y = yBaseGrid + PARAM.CutStep*ny.*sum(repmat(perturb,1,numel(x)).*PPP(2:PARAM.dealiasing,:))';
        elseif PARAM.dropFrame==1
            x = xBaseGrid + PARAM.CutStep*nx.*sum(repmat(perturb,1,numel(x)).*PPP(3:PARAM.dealiasing,:))';
            y = yBaseGrid + PARAM.CutStep*ny.*sum(repmat(perturb,1,numel(x)).*PPP(3:PARAM.dealiasing,:))';
        end
        
        %dealiasing
        [x,y] = dealiasingGridXY(x,y,PARAM);
        
        %update ghost point
        if PARAM.dropFrame==0
            fVolume = @(unk) ModifyVolumeSpectralXYdns(x,y,nx,ny,unk,V0,PARAM);
            options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
            move = fsolve(fVolume,0,options);
            xBaseGrid = x + nx*move;
            yBaseGrid = y + ny*move;
        elseif PARAM.dropFrame==0
            fVolume = @(unk) ModifyVolumeCmSpectralXY(x,y,nx,ny,unk,V0,xcm,PARAM);
            options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
            move = fsolve(fVolume,[0 0],options);
            xBaseGrid = x + nx*move(1).*PPP(1,:)' + nx*move(2).*PPP(2,:)';
            yBaseGrid = y + ny*move(1).*PPP(1,:)' + ny*move(2).*PPP(2,:)';
        end
        
        %dealiasing
        [xBaseGrid,yBaseGrid] = dealiasingGridXY(xBaseGrid,yBaseGrid,PARAM);
        
        %update solution and base modes
        sol(1:end-1) = 0;
        
        %compute the modes from the grid points
        [xBase,yBase] = fromGridToModes(xBaseGrid,yBaseGrid,PARAM);
        
        %dealiasing
        xBase = xBase(1:PARAM.dealiasing);
        yBase = yBase(1:PARAM.dealiasing);
        
        %nonlinear equation
        fNonlinear = @(unk) NormalVelocitySpectralVolumeXYmodesCont(unk,xBase,yBase,V0,xcm,PARAM,sol(end));
        [~,nx,ny,xBaseGrid,yBaseGrid] = fNonlinear(sol(1:end-1));
        
        solPrev = sol;
        sol = sol+dir*delta; % new prediction of solution
        
    else
        
        solPrev = sol;
        sol = sol+dir*delta; % new prediction of solution
        
    end
    
    %remesh
    %if remesh==1 && sum(0:howOftenRemesh:2000==countCont) && minEig>isSing
    
    %compute modes
    [xMode,yMode] = fromGridToModes(x,y,PARAM);
    xyMode(1:2:end-1) = xMode;
    xyMode(2:2:end) = yMode;
    
    if PARAM.remesh==1 && checkMetricCurvilinear(0,xyMode,PARAM) && minEig>isSing && countCont>1
        
            disp('Remesh')
            xBase(PARAM.dealiasing+1:numel(x)) = 0;
            yBase(PARAM.dealiasing+1:numel(x)) = 0;
            yyy = zeros(2*numel(xBase),1);
            yyy(1:2:end-1) = xBase; yyy(2:2:end) = yBase;
            [yyy,halt] = remeshDropSpectral(yyy,PARAM);
            if halt==1
                warning('Remesh has failed');
            elseif halt==0
                xBase = yyy(1:2:end-1); yBase = yyy(2:2:end);
                xBase = xBase(1:PARAM.dealiasing);
                yBase = yBase(1:PARAM.dealiasing);

                %compensate previous augmentation
                sol(1:end-1) = sol(1:end-1)-dir(1:end-1)*delta;    sol(end) = manyCa(countCont-1)-delta*sign(dir(end));

                %nonlinear equation
                fNonlinear = @(unk) NormalVelocitySpectralVolumeXYmodesCont(unk,xBase,yBase,V0,xcm,PARAM,sol(end));
                [~,nx,ny,xBaseGrid,yBaseGrid] = fNonlinear(sol(1:end-1));

                solPrev = sol;
                if PARAM.dropFrame==0
                    dir = [zeros(PARAM.dealiasing-1,1); sign(dir(end))];
                elseif PARAM.dropFrame==1
                    dir = [zeros(PARAM.dealiasing-2,1); sign(dir(end))];
                end
                sol = sol+dir*delta;
            end
            
    end

    display(['Newton loop ' num2str(countCont) ', Ca=' num2str(sol(end)) ' lambda=' num2str(PARAM.visc)])

    % Newton iterations
    quit=0;count=0;
    while ~quit

        %present solution
        perturb = sol(1:end-1); PARAM.Ca = sol(end);
        
        %nonlinear equation
        fNonlinear = @(unk) NormalVelocitySpectralVolumeXYmodesCont(unk,xBase,yBase,V0,xcm,PARAM,PARAM.Ca);
        fNonlinearCa = @(unk) NormalVelocitySpectralVolumeXYmodesCont(perturb,xBase,yBase,V0,xcm,PARAM,unk);

        %NON LINEAR SOLUTION
        %normal velocity
        [u,~,~,x,y] = fNonlinear(perturb);
        R(1:end-1) = u;
        %continuation
        R(end) = dir'*(sol-solPrev)-delta;
        
        if PARAM.computeOnlyOneJacobian==1
            disp('Compute only on Jacobian')
        end
    
        %Jacobians
        if PARAM.cpu==1
            %Jacobian
            if (count==0 && sum(1:PARAM.computeJacobianHowOften:1e4==countCont) && PARAM.computeOnlyOneJacobian==1) || count>PARAM.computeJacobianNow
                J = JacobianHandle(fNonlinear,perturb,PARAM.dh);
            elseif PARAM.computeOnlyOneJacobian==0
                J = JacobianHandle(fNonlinear,perturb,PARAM.dh);
            end
            %J = JacobianHandle(fNonlinear,perturb,PARAM.dh);
        elseif PARAM.cpu>1
            J = JacobianHandleParallel(fNonlinear,perturb,PARAM.dh);
        end
        Jca = JacobianHandle(fNonlinearCa,PARAM.Ca,PARAM.dh);
        
        %build Jacobian
        JJJ = [J Jca; dir'];

        %compute elongation
        D = max(x);
        %AAA = x(1)-x(end);
        %BBB = 2*y(round(numel(y)/2));
        %D = (AAA-BBB)/(AAA+BBB);

        %new guess
        step = JJJ\R;
        sol = sol - step*PARAM.CutStep;

        % convergence test
        res = max(abs(R));
        
        figure(1)
        plot([x; flip(x)],[y; -flip(y)],[aIN; flip(aIN)],[bIN; -flip(bIN)],'--');
        grid on; axis equal; axis([-4 4 -2 2]);%axis([1.9 2.1 -0.05 0.05]);
        disp(['Iteration ' num2str(count+1) ' res=' num2str(res) '  D=' num2str(D)]);
        legend('Solution','Initial guess','Location','Best')
        title('Droplet shape')
        if count>PARAM.stop||res>1e1||isnan(res)==1; disp('no convergence'); break; end
        if res < PARAM.ResBreak; quit=1; disp('converged'); continue; end

        %count newton iteration
        count = count+1;

    end
    
    %compute eigenvalues
    [eigenMode,eigenval] = eig(J);
    eigenval = diag(eigenval);
    unstEigen(countCont) = sum(real(eigenval)>0)-1;
    colEigen = sum(real(eigenval)>0);
    if colEigen==0
        colEigen = 5;
    end
    colHere = [col(colEigen) 'o'];
    
    figure(2)
    if countCont>2
    hold on
    end
    plot(sol(end),D,colHere)
    hold off
    grid on
    title(['Bifurcation diagram \lambda=' num2str(PARAM.visc)])
    xlabel('Ca')
    ylabel('max(x)')
    drawnow
    
    %break if there is no convergence
    if PARAM.legendre==1
        onAxis = sum((y<0)>0);
    elseif PARAM.legendre==0||PARAM.legendre==2
        onAxis = sum((y(2:end-1)<-1e-4)>0);
    end
    
    if res>1e1 || isnan(res)==1 || count>PARAM.stop || onAxis
        disp('Newton method did not converged or interface has crossed the axis')
        break
    end

    manyD(countCont) = D;
    manyCa(countCont) = PARAM.Ca;
    manyA(:,countCont) = x;
    manyB(:,countCont) = y;
    manyPerturb(:,countCont) = perturb;
    manyXbase(:,countCont) = xBase;
    manyYbase(:,countCont) = yBase;
    %storeJacobian{countCont} = J;

    cd(PARAM.res)
    save(PARAM.filename)
    cd(here)

    % New direction
    if PARAM.dropFrame==0
        newDirextionVector = [zeros(PARAM.dealiasing-1,1); 1];
    elseif PARAM.dropFrame==1
        newDirextionVector = [zeros(PARAM.dealiasing-2,1); 1];
    end
    dir=JJJ\newDirextionVector; % nouvelle direction
    dir=dir/norm(dir); % normalization
    
    %figure out if the Jacobian is singular
    [minEig,indMin] = min(abs(eigenval));
    if minEig<isSing && beforeSing==0
        disp('The Jacobian is singular')
        beforeSing = beforeSing+1;
        
        if switchBranch==1
            
            disp('Switch branch')
           
            %find neutral mode
            neutralMode = eigenMode(:,indMin)/max(eigenMode(:,indMin));
            
            %add the neutral mode to the solution
            sol(1:end-1) = sol(1:end-1) + deltaSwitch*neutralMode;
            dir(end) = dir(end)*deltaDir;
            dir=dir/norm(dir); % normalization
            %dir = dir*deltaDir;
            
            if checkOnlyOnePitchfork==1
                switchBranch = 0;
            end
            
        end
        
    end
    
    if beforeSing>1
        beforeSing = 0;
    end

    countCont = countCont+1;

end

cd(PARAM.res)
disp('Save results')
save(PARAM.filename)
cd(here)