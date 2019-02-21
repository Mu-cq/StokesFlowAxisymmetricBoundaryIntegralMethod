%compute interface shape for droplet in extensional flow usinf Newton
%method

function minimalSeedOneDropBEM(Tsave,maxDT,initialDT,PARAM,V0,A0perturb)

error('Still in development')

printToScreenOneDrop(PARAM.Ca,PARAM.Bond,PARAM.visc,PARAM.dropFrame,PARAM.algorithm);

%current directory
here = pwd;

%initial shape
[initialXY,PARAM] = initialConditionDropBEM(PARAM);
xBase = initialXY(1:2:end-1);
yBase = initialXY(2:2:end);
areaBaseState = surf_gauss_vect(xBase',yBase');

%plot initial guess
figure
plot(xBase,yBase,'k')
hold on
axis equal
plot(xBase,-yBase,'k')
xlabel('z')
ylabel('r')
title('droplet shape')
grid on
drawnow

%solve non linear optimization problem subject to nonlinear constraint
A = []; b = []; Aeq = []; beq = []; lb = -10*ones(1,numel(initialXY)); ub = 10*ones(1,numel(initialXY));
%lb(2:2:end) = -0.01;

%function handles
fToMinimize = @(xyUNK) runTimeSteppingMinimalSeedBEM(xyUNK,Tsave,maxDT,initialDT,PARAM,V0);  % nonlinear equation
fConstr = @(xyUNK) conserveSurfaceAreaBEM(xyUNK,V0,A0perturb+areaBaseState);                   % nonlinear constraints

%check initial shape
if min(yBase)<0
    warning('Negative value of y coordinate is unphysical')
end

%options
%options = optimoptions('fmincon','Display','iter','Algorithm','sqp','TolCon',1e-06,'TolX',1e-06,'TolFun',1e-06,'MaxFunEvals',1e4);

%minimize
XY = fmincon(fToMinimize,initialXY,A,b,Aeq,beq,lb,ub,fConstr);

%get final shape
x = XY(1:2:end-1);  y = XY(2:2:end);

%plot final shape
%figure
plot(x,y,'r')
hold on
axis equal
plot(x,-y,'r')
xlabel('z')
ylabel('r')
title('droplet shape')
grid on
drawnow

% if PARAM.plotCurv==1
%     %plot curvature
%     figure(5)
%     plot(x,K1,'-')
%     grid on
%     hold on
%     plot(x,K2,'-')
%     plot(x,K1+K2,'-k')
%     title(['Curvature Ca=' num2str(PARAM.Ca)])
%     legend('K_1','K_2','K_1+K_2','Location','Best')
%     xlabel('x')
%     ylabel('K')
% end

cd(PARAM.res)
save(PARAM.filename)
cd(here)