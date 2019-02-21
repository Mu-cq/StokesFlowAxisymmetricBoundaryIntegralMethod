%print to screen information when starting the simulation

function printToScreenTimeStepping(ODE,dt,tol)

if ODE==1
    ODEout = 'ODE45';
elseif ODE==2
    ODEout = 'ODE23t';
elseif ODE==3
    ODEout = 'ODE23s';
elseif ODE==0
    ODEout = 'RK2';
elseif ODE==4
    ODEout = 'ODE15s';
end

%printf(' \n',)

display(['Time stepping is ' ODEout ' with dt=' num2str(dt) ' and tol=' num2str(tol)])






