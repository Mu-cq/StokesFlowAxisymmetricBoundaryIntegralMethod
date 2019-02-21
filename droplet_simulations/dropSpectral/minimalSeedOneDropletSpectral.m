%compute interface shape for droplet in extensional flow usinf Newton
%method

function minimalSeedOneDropletSpectral(Tsave,maxDT,initialDT,PARAM)

error('In developement')

tic

%current directory
here = pwd;

%print to screen
PARAM.Tend = Tsave(end);
printToScreen(PARAM);

%final volume that I want
V0 = PARAM.V0;

%initial (base) shape
if PARAM.BC==1
    
    if PARAM.uploadShape==2 || PARAM.uploadShape==5
        %upload shape
        [~,xBase,yBase] = initialConditionDrop(PARAM);

        %upload shape from Newton in order to compute the surface area base
        %state
        [xBaseTemp,yBaseTemp] = uploadFromNewton(PARAM);
        areaBaseState = surfaceCurvilinearAxisSpectral(xBaseTemp,yBaseTemp,PARAM);

    else
        error('You need to start from a base state obtained with Newton method or Minimal Seed')
    end

elseif PARAM.BC==2

    %the base state is known analitically
    r0 = nthroot(3/4/pi*PARAM.V0,3);
    areaBaseState = 4*pi*r0^2;
    
    if PARAM.uploadShape==0 && PARAM.D~=0
        error('Base state for rising droplet is a sphere')
    end
    
    %upload shape
    [~,xBase,yBase] = initialConditionDrop(PARAM);
    
end
%plot initial shape
% figure
% plot([xBase; flip(xBase)],[yBase; -flip(yBase)],'k')
% hold on
% axis equal
% axis ([-2 2 -2 2])
% %plot(xBase,-yBase,'k')
% xlabel('z')
% ylabel('r')
% title('droplet shape')
% grid on
% drawnow

if PARAM.minimalSeedOptim==1    % use fMinCon
    
    error('Very often it does not converge')
    
    %initialize
    perturb = zeros(PARAM.dealiasing-1,1);

    %solve non linear optimization problem subject to nonlinear constraint
    A = []; b = []; Aeq = []; beq = []; lb = -10*ones(1,2*numel(perturb)); ub = 10*ones(1,2*numel(perturb));

    %function handles
    fToMinimize = @(perturb) runTimeSteppingMinimalSeed(perturb,xBase,yBase,Tsave,maxDT,initialDT,PARAM,V0);  % nonlinear equation
    fConstr = @(perturb) conserveSurfaceAreaSpectral(perturb,xBase,yBase,V0,PARAM.A0perturb+areaBaseState,PARAM);                   % nonlinear constraints
    [~,ny] = normalVectorSpectral(xBase,yBase,PARAM);

    %check initial shape
    if min(ny)<-1e-5
        warning('Normal vector of the initial shape is not accurate')
    end

    %options
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp','TolCon',1e-06,'TolX',1e-06,'TolFun',1e-06,'MaxFunEvals',1e4);

    %minimize
    %[pertub,fval,exitflag,output] = fmincon(fToMinimize,XY,A,b,Aeq,beq,lb,ub,fConstr,options);
    perturb = fmincon(fToMinimize,perturb,A,b,Aeq,beq,lb,ub,fConstr,options);
    
    %get final shape
    [x,y] = getVolumeConservingShapeSpectral(perturb,xBase,yBase,V0,PARAM);

elseif PARAM.minimalSeedOptim==2    % use gradient descent, in house
    
    %initialize
    perturb = zeros(PARAM.dealiasing-2,1);
    
    %objective function
    fToMinimize = @(perturb) runTimeSteppingMinimalSeedConserveArea(perturb,xBase,yBase,Tsave,maxDT,initialDT,PARAM,V0,PARAM.A0perturb+areaBaseState);  % nonlinear equation
    
    %start optimization
    OutFun0 = zeros(PARAM.iterMinimalSeed,1);
    T = cell(PARAM.iterMinimalSeed,1);
    Y = cell(PARAM.iterMinimalSeed,1);
    for i = 1:PARAM.iterMinimalSeed
        
       OutFun = zeros(numel(perturb),1);
       [OutFun0(i),T{i},Y{i}] = fToMinimize(perturb);
       %compute objective function for every varaition of DOF
       for k = 1:numel(perturb)
           
          perturbHere = perturb;
          perturbHere(k) = perturbHere(k)+PARAM.dh;
          OutFun(k) = fToMinimize(perturbHere);
           
       end
       
       %compute Jacobian
       Jacobian = (OutFun-OutFun0(i))/PARAM.dh;
       
       %move toward descent
       perturb = perturb - PARAM.stepOptim*Jacobian;
       
       %display status of the iteration
       display(['Gradient descent: iteration=' num2str(i) ' of ' num2str(PARAM.iterMinimalSeed) ', objective function F=' num2str(OutFun0(i)+areaBaseState)])
        
    end
    
    %plot final shape
    [x,y] = getVolumeAndAreaConservingShapeSpectral(perturb,xBase,yBase,V0,PARAM.A0perturb+areaBaseState,PARAM);
    
end

%plot final shape
%figure
% plot([x; flip(x)],[y; -flip(y)],'r')
% %plot(x,-y,'r')
% xlabel('z')
% ylabel('r')
% title('droplet shape')
% legend('Initial shape','Final shape','Location','Best')
% grid on
% drawnow

simulationTime = toc;

cd(PARAM.res)
display('Save results')
save(PARAM.filename)
cd(here)
