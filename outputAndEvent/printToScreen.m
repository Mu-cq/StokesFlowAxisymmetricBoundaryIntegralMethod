%print to screen information when starting the simulation

function printToScreen(PARAM)

%boundary condition
if PARAM.BC==1
    type = 'drop in linear extensional flow';
elseif PARAM.BC==2
    type = 'rising droplet';
elseif PARAM.BC==3
    type = 'drop in nonlinear extensional flow';
end

%frame of reference
if PARAM.dropFrame==0 || PARAM.dropFrame==2
    frame = 'lab';
elseif PARAM.dropFrame==1
    frame = 'co-moving';
end

% time integration
ode = '';
if PARAM.ODE==1
    ode = 'ODE45';
elseif PARAM.ODE==6
    ode = 'ODE23t';
elseif PARAM.ODE==2
    ode = 'myRK2';
elseif PARAM.ODE==9
    ode = 'myRK4';
end

%spetrcal basis
spectral = '';
if PARAM.legendre==1
    spectral = 'Legendre';
elseif PARAM.legendre==0
    spectral = 'Chebyshev';
elseif PARAM.legendre==2
    spectral = 'Legendre-Lobatto';
end

%algorithm
if PARAM.algorithm==1

    display(['One ' type ' starts in the ' frame ' frame with ' num2str(PARAM.dealiasing) ' ' spectral ' modes, time stepper is ' ode ' Ca=' num2str(PARAM.Ca) ' visc=' num2str(PARAM.visc) ])

elseif PARAM.algorithm==2
    
    display(['Edge tracking one ' type ', loop number ' num2str(PARAM.loopNum) ' starts in the ' frame ' frame with ' num2str(PARAM.dealiasing) ' ' spectral ' modes, time stepper is ' ode ' Ca=' num2str(PARAM.Ca) ' visc=' num2str(PARAM.visc) ])
    
elseif PARAM.algorithm==3
    
    display(['Newton method one ' type ' starts in the ' frame ' frame with ' num2str(PARAM.dealiasing) ' ' spectral ' modes, Ca=' num2str(PARAM.Ca) ' visc=' num2str(PARAM.visc) ])
    
elseif PARAM.algorithm==4
    
    display(['Continuation one ' type ' starts in the ' frame ' frame with ' num2str(PARAM.dealiasing) ' ' spectral ' modes, delta=' num2str(PARAM.delta) ', CaUp=' num2str(PARAM.CaBreakUp) ', CaDown=' num2str(PARAM.CaBreakDown) ', visc=' num2str(PARAM.visc) ])
    
elseif PARAM.algorithm==5
    
    display(['Stability analysis in the ' frame ' frame for one ' type ' Ca=' num2str(PARAM.Ca)])
    
elseif PARAM.algorithm==6
    
    display(['Minimal seed for break-up in the ' frame ' frame for one ' type ' Ca=' num2str(PARAM.Ca)])
    display(['Energy of the initial perturbation is A0=' num2str(PARAM.A0perturb) ', and time horizon is T=' num2str(PARAM.Tend)])
    
else
    
    error('Invalid algorithm selected')

end