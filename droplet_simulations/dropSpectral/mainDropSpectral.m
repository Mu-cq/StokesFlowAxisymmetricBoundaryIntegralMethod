%nonlinear BEM based solver for one droplet motion
%spectral with Chebyshev or Legendre

clear variables
close all

%save data
here = pwd;
results = '~/Documents/MATLAB/droplet_simulations/results';

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.algorithm = 3;               % 1 is time marching DNS, 2 is edge tracking, 3 is newton method, 4 is continuation, 5 is stability analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% PHYSICAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.Ca = 0.1;                                       % capillary number
PARAM.visc = 1;     PARAM.Qdeflation = 0;           % viscosity ratio and deflation
V0 = 4/3*pi;                                        % drop volume
PARAM.placeShapeXCM = 0;                            % when doing edge tracking, replace always droplet in x0=0 at t0
PARAM.CaNL = 0;                                     % non linear extensional flow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose where to get the shape
PARAM.uploadShape = 0;                                       % 0 the initial shape is analytical, 1 is from continuation, 2 the initial shape is from newton Method,  3 is from edge tracking, 4 is from DNS, 5 is from Minimal Seed
PARAM.D = 0.1;                                              % deformation parameter of the initial condition

% in case the initial shape is analytical, choose how to draw it
PARAM.shapeEllipse = 0;                                      % 0 initial shape is sphere plus P2 Legendre, 1 initial shape is an ellipse, 2 is pertubed with 2 legendre polynomials P2 P4, 3 is asymmetric shape by P3 Legendre, 4 is P2, P3
PARAM.f2 = 0.8;                                              % if I use first two leegndre polynomials

% if the shape is uploaded from previous simulation
PARAM.elemUpload = 50;                                       % number of elemnt of uploaded data
PARAM.BCupload = 2;                                          % BC of the uploaded shape
PARAM.legendreUpload = 1;                                    % spectral basis of the previous simulation
PARAM.overlapMode = 1;  PARAM.whichmode = 2;                    % 1 perturb with legendre polynomial, 2 perturb shape with eigenmode

% if it start from continuation, newton or edge tracking
PARAM.CaUpUpload = 1; PARAM.CaDownUpload = -1;              % for uploading the file
PARAM.Dupload = 0;                                           % choose defromation parameter on bifurcation diagram
PARAM.CaUpload = 0.25;                                     % choose Ca on bifurcation diagram
PARAM.viscUpload = PARAM.visc;                                  % choose viscosity ratio for upload

% if the shape is uploaded from edge tracking or DNS
PARAM.edgeLoopUpload = 1000;                                   % loops of the edge tracking simulation
PARAM.deltaEdgeUpload = 1e-4;                                % delta edge of the uloaded data
PARAM.dtUpload = 2e-4;                                       % time step of the uploaded simulation
PARAM.volCorrUpload = 1;                                     % if the simulation was using volume correction
PARAM.whichLoopUpload = 765;                                  % which loop of the edge tracking has to be loaded
PARAM.Tupload = 26.85;                                         % choose at which time upload shape

% if shape is uploaded from DNS
PARAM.TendUpload = 100;                                         % choose at which time upload shape

% if shape is uploaded from Minimal Seed
PARAM.A0upload = 0.2;                                         % energy amplification of previous optimization
PARAM.ThorizonUpload = 10;                                   % time horizon of previous optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.dropFrame = 0;                                % choose frame of reference: 0 is lab frame, 1 is drop frame
PARAM.Unormal = 1;                                  % 0 advect nodes with lagrangian velocity 1 with normal velocity, 3 use maesh stab, 4 radial velocity
PARAM.volume = 0;   PARAM.volTol = 1e-14;           % impose volume conservation with first mode
PARAM.BC = 1;                                       % BC: 1 is linear extensional flow, 2 is rising droplet, 3 is nonlinear extensional flow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REMESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.remesh = 1;           % if 1 remesh is activated with newton method (spectral accuracy), 2 is with spline (lose spectral accuracy)
PARAM.remeshUniform = 0;    % impose uniform distribution of the mesh
PARAM.chooseMapping = 1;  PARAM.howStrong = 0.85; % 1 is normal mapping (like spectral parametrization), 2 is clustered toward the middle (how strong set how stronly to cluster)
PARAM.tolRemesh = 1e-5;     % tolerance for remeshing
PARAM.firstTimeRemesh = 0;  % don't perform remesh before this instant
PARAM.normRemesh = 1;       % norm on which I base the remeshing criteria: 1 is on the first derivative, 2 is on  the second derivative, 3 is on the sum of the two
PARAM.remeshSolve = 1;      % 1 newton method or 2 minimization fot remeshing
PARAM.MaxIter = 50;         % maximum iteration of "fsolve" or "fmincon"
PARAM.remeshStart = 0;      % remesh when drawning initial shape
PARAM.checkRemeshFail = 1;  % check if newton method for remeshing fails

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%% SPACE DISCRETIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.legendre = 1;                                         % choose 0 Chebyshev, 1 Legendre, 2 is Legendre-Lobatto
PARAM.dealiasing = 50;  ratioGrid = 2;                     % dealiasing
PARAM.tresholdCoeffs = 1e-2;     PARAM.Rescale = 1;         % option for accuracy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% TIME DISCRETIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.ODE = 2;          % 1 is ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
Tstart = 0;             % beginning of simulation
Tend = 15;             % end of simulation
maxDT = 2e-4;           % set max time step
initialDT = maxDT;      % set initial time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% EDGE TRACKING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%

PARAM.edgeLoop = 100;                          % loop of edge tracking
PARAM.edgeStartingLoop = 1;                     % starting loop of edge tracking
PARAM.bisection = 2;                            % 1 average points resizing for volume conservation, 2 as before but resize also to bisect the area value
PARAM.deltaEdge = 1e-2;                         % next loop starts when a norm is smaller that this
PARAM.normEdge = 2;                             % on which norm I do the check: 1 is L-2 norm, 2 id L-infinity norm, 3 is excess area norm
PARAM.Dedge = [0.8 1];                       % deformation parameter first two shapes to check
PARAM.overlapModeEdge = [1 1];                  % kind of mode for perturbing first two shapes to check
PARAM.whichModeEdge = [2 2];                    % exactly which mode
PARAM.bruteForceBisection = 0;                  % ingnore bisection criterion, just bisect between two previous shapes
PARAM.deltaModify = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% SAVING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.saveData = 1;                                 % choose if to save
PARAM.saveDataIte = 0;                              % choose if to save every iteration
SaveHowMany = 150;                                  % output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERGENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%residuals options
PARAM.checkRes = 1;                 % check residuals or not
PARAM.converge = 1e-6;              % when resilduals are smaller, convergence
PARAM.convergeShape = -0.01;         % convergence based on shape
PARAM.convergeShapeEdge = 0.02;     % convergence based on shape

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEWTON METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.plotRes = 1;                      % plot residuals
PARAM.ResBreak = 1e-10;                 % break when residuals are smaller
PARAM.stop = 200;                       % stop after maximum iteration
PARAM.CutStep = 1;                      % advance only partially
PARAM.plotCurv = 1;                     % plot curvature
PARAM.dh = 1e-5;                        % finite diffrence for Jacobian computation
PARAM.computeOnlyOneJacobian = 0;       % compute the Jacobian only at the first iteration
PARAM.call_fsolve = 0;                  % use matlab fsolve function
PARAM.parallel_jacobian_on = 1;         % use parallel toolbox for the computation of the Jacobian
PARAM.cpu = [];                         % number of cpu for parallel jacobian, is [] take workers by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONTINUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.delta = 0.01;                                                             % continaution method step
PARAM.CaBreakUp = 1;                                                         % break when Ca is larger
PARAM.CaBreakDown = 0.07;                                                        % break when Ca is smaller
PARAM.changeSystemOfCoord = 1;    PARAM.howOften = 1;                           % change system of coordinates every continuation iteration
PARAM.switchBranch = 1;   PARAM.deltaSwitch  = 0.1;    PARAM.deltaDir = 1;      % for branch switching
PARAM.isSing = 1e-3;      PARAM.checkOnlyOnePitchfork = 1;
PARAM.computeJacobianHowOften = 1;                                             % how often the Jacobian is computed
PARAM.computeJacobianNow = PARAM.stop/2;                                        % if Newton does not converge, compute the Jacobian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MINIMAL SEED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.A0perturb = 0.09;                    % enery of the perturbation
PARAM.minimalSeedOptim = 2;             % 1 is fMinCon (not working well, probably the shape is changing too much compared to the base), 2 is gradient descent (in house)
PARAM.iterMinimalSeed = 20;             % number of iteration for the optimization algorithm 
PARAM.stepOptim = 0.1;                 % step for the Jacobian during optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of grid points
PARAM.n = round(ratioGrid*PARAM.dealiasing)+1;
PARAM.V0 = V0;
PARAM.Tend = Tend;
if Tstart~=PARAM.Tupload && PARAM.uploadShape==4
    error('Upload and start time have to be the same');
end

%filename
PARAM.filename = chooseFilename(PARAM,maxDT,Tend,PARAM.A0perturb);
PARAM.filenameIte = PARAM.filename;
PARAM.res = results;
PARAM.this_folder = here;

%choose spectral mapping
PARAM.remeshMapping = chooseSpectralMapping(PARAM.chooseMapping,PARAM.howStrong);

%compute integration weight and differentiation matrices
if PARAM.legendre==1
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreIntDiff(0,1,PARAM.n+1);
    [PARAM.tcheb,PARAM.D1cheb,PARAM.D2cheb,PARAM.WGcheb] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==0
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.TTT] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==2
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreLobattoIntDiff(0,1,PARAM.n+1);
end

if PARAM.parallel_jacobian_on==1 && PARAM.call_fsolve==1
    disp('Compute Jacobian in parallel')
    p = gcp('nocreate');
    if isempty(p)==1
        if isempty(PARAM.cpu)==0
            parpool(PARAM.cpu)
        else
            parpool()
        end
    end
end
  
%run computation
if PARAM.algorithm==1
    %time stepping
    [T,Y] = runTimeStepping(Tsave,maxDT,initialDT,PARAM,V0);
elseif PARAM.algorithm==2
    %edge tracking
    [Tedge,Yedge] = edgeTracking(Tsave,maxDT,initialDT,PARAM,V0);
elseif PARAM.algorithm==3
    %newton method
    newtonMethodSpectralXYmodesFunction(PARAM);
elseif PARAM.algorithm==4
    %newton method
    newtonMethodSpectralXYmodesContFunction(PARAM);
elseif PARAM.algorithm==5 && (PARAM.BC==1 || PARAM.BC==3)
    %stability analyisis droplet in extensioanl flow
    stabilitySpectralXYmodesFunction(PARAM)
elseif PARAM.algorithm==6
    %minimal seed for droplet break-up
    %error('Problem during the convergence')
    minimalSeedOneDropletSpectral(Tsave,maxDT,initialDT,PARAM);
end

%take simulation time
simulationTime = toc;

%save results
if PARAM.saveData==1 && PARAM.algorithm==1
    cd(results)
    disp('Save data')
    save(PARAM.filename)
    cd(here)
end

%the end
disp('The End')















