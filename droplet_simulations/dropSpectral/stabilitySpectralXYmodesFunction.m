%stability of droplet interface by numerical perturbing the Jacobian

function stabilitySpectralXYmodesFunction(PARAM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotModes = 2;                      % 0 does not plot eigenmodes, 1 plot all eigemodes, 2 plot unstable eigemodes
overlap = 1;    ampli = -0.2;        % overlap unsatbale eigenmodes to base shape with arbitrary amplitude
eigIsPos = 1e-5;                    % consider eigenvalue as neutral

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%upload shape
PARAMold = PARAM;
[~,~,~,newCa,xBase,yBase,perturb,PARAM] = initialConditionDrop(PARAM);
if PARAMold.uploadShape==1
    PARAMold.Ca = newCa;
    PARAM.Ca = newCa;
end

if PARAMold.dropFrame==1
    if PARAM.dropFrame==0
        error('Previous simulation was in a different frame of reference')
    end
    PARAM.dropFrame = 1;
elseif PARAMold.dropFrame==0
    if PARAM.dropFrame==1
        error('Previous simulation was in a different frame of reference')
    end
    PARAM.dropFrame = 0;
end

%ofetn used variables
V0 = PARAM.V0;
[x,y] = fromModesToGrid(xBase,xBase,PARAM);
xcm = CenterMassCurvAxisSpectral(x,y,PARAM);

%print output
printToScreen(PARAMold);

%nonlinear equation
fNonlinear = @(unk) NormalVelocitySpectralVolumeXY2modes(unk,xBase,yBase,V0,xcm,PARAM);

%check that this is actualy a solution
[u,nx,ny,x,y,uFirst,K1,K2] = fNonlinear(perturb);

%plot shape which is a solution
figure(1)
%subplot(2,2,1)
plot([x; flip(x)],[y; -flip(y)],'k')
axis([-3 3 -3 3])
xlabel('$z$','interpreter','latex')
ylabel('$r$','interpreter','latex')
grid on
axis equal
title(['Shape Ca=' num2str(PARAM.Ca) ' \lambda=' num2str(PARAM.visc)])

%plot modes coefficients
figure(5)
%subplot(2,2,2)
loglog(abs(xBase(2:2:end)))
grid on
hold on
loglog(abs(yBase(1:2:end-1)))
title('Modes coefficients')
legend('x','y','Location','Best')
xlabel('x')
ylabel('C_n')

%plot curvature
figure(2)
%subplot(2,2,2)
plot(x,K1,'-')
grid on
hold on
plot(x,K2,'-')
plot(x,K1+K2,'-k')
title('Curvature of the base state')
legend('K_1','K_2','K_1+K_2','Location','Best')
xlabel('x')
ylabel('K')

if max(abs([u; uFirst']))>PARAM.ResBreak
    warning('Current shape it is not a solution because the residuals are too high')
    figure(10)
    loglog(abs(u(2:2:end)))
    grid on
    title('Velocity coefficients')
    xlabel('x')
    ylabel('U_n')
end
if abs((VolumeCurvilinearAxisSpectral(x,y,PARAM)-V0)/V0)>1e-12
    warning('Current shape is not of the right volume')
end
if ny<0
    warning('Current shape it is not accurate because ny is negative')
end

%compute numerically the Jacobian
J = JacobianHandle(fNonlinear,perturb,PARAM.dh);

%compute eigenvalues and eigenfunctions
[eigenmode,eigenval] = eig(J);
eigenval  = diag(eigenval);

%order from small to big
[eigenvalR,ind] = sort(real(eigenval),'descend');
eigenval = eigenvalR + 1i*imag(eigenval(ind));
eigenmode = eigenmode(:,ind);

%plot eigenvalues and eigenfunction
modePhysical = zeros(numel(PARAM.t),numel(eigenval));
for i = 1:numel(eigenval)
    
    figure(3)
    %subplot(2,2,3)
    plot(real(eigenval(i)),imag(eigenval(i)),'.','MarkerSize',30)
    hold on
    xlabel('$\sigma_r$','interpreter','latex')
    ylabel('$\sigma_i$','interpreter','latex')
    grid on
    title('Eigenvalues')
    
    if plotModes>=1
    
        %compute physical eigenmode
        mode = eigenmode(:,i);
        %compute current grid points
        if PARAM.legendre==1
            PPP = PARAM.PPP;
            %modePhysical(:,i) = sum(repmat(mode,1,numel(x)).*PPP(2:PARAM.dealiasing,:));
        elseif PARAM.legendre==0
            PPP = PARAM.TTT;
            %modePhysical(:,i) = sum(repmat(mode,1,numel(x)).*TTT(2:PARAM.dealiasing,:));
        end
        
        if PARAM.dropFrame==0
                modePhysical(:,i) = sum(repmat(mode,1,numel(x)).*PPP(2:PARAM.dealiasing,:));
        elseif PARAM.dropFrame==1 || PARAM.dropFrame==2
                modePhysical(:,i) = sum(repmat(mode,1,numel(x)).*PPP(3:PARAM.dealiasing,:));
        end
    
        if plotModes==1
            plotYes = 1;
        elseif plotModes==2 && real(eigenval(i))>eigIsPos
            plotYes = 1;
        else
            plotYes = 0;
        end
        
        if plotYes==1
        
            figure(4)
            %subplot(2,2,4)
            plot(PARAM.t,ampli*modePhysical(:,i)/modePhysical(1,i))
            %plot(PARAM.t,modePhysical(:,i)/modePhysical(1,i))
            hold on
            grid on
            title('Eigenmodes')
            xlabel('$s$','interpreter','latex')
            ylabel('$\hat{\eta}$','interpreter','latex')
            
        end
    
    end
    
end

%number of positive eigenvalues
positive = sum(real(eigenval)>eigIsPos);
if positive==1
    disp(['There is ' num2str(positive) ' positive eigenvalue'])
else
    disp(['There are ' num2str(positive) ' positive eigenvalues'])
end

if overlap==1
   
    for k = 1:positive
        
        %perturb respect to base shape
        overlapMode = ampli*modePhysical(:,k)/max(modePhysical(:,k));
        xEigen = x + nx.*overlapMode;
        yEigen = y + ny.*overlapMode;
    
        figure(1)
        %subplot(2,2,1)
        hold on
        plot([xEigen; flip(xEigen)],[yEigen; -flip(yEigen)])
        axis([-3 3 -3 3])
        
        %save(['xEigenCa=' num2str(PARAM.Ca) '.mat'],'xEigen')
        %save(['yEigenCa=' num2str(PARAM.Ca) '.mat'],'yEigen')
    
    end
    
end