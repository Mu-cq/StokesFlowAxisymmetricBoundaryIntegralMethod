%Runge-Kutta order 2

function [T,X,Y,U,V,laplace,stokes] = RK2panels(f,TTT,initial,DT,initialDT,remesh,event)

%take variables step
if initialDT==DT
    ttt = TTT(1):DT:TTT(end);
else
    ttt = [TTT(1):initialDT:TTT(20)-initialDT TTT(20):DT:TTT(end)];
end
manyDT = diff(ttt);
loop = numel(ttt);

%check dt and saving time
if sum(sum(((repmat(ttt,numel(TTT),1)>repmat(TTT',1,numel(ttt))-1e-6) + (repmat(ttt,numel(TTT),1)<repmat(TTT',1,numel(ttt))+1e-6))==2))<numel(TTT)
    warning('Saving frequency is not a correct multiple of the time step, saving steps can be lost')
end

%initialize
X = cell(numel(TTT),1);
U = cell(numel(TTT),1);
Y = cell(numel(TTT),1);
V = cell(numel(TTT),1);
laplace = cell(numel(TTT),1);
stokes = cell(numel(TTT),1);
T = zeros(numel(TTT),1);

%initial condition
x2 = initial{1};
y2 = initial{2};
tParametric = initial{3};
PARAM_2 = initial{4};

count = 1;
for i = 1:loop
    
    %current time
    t =ttt(i);
    
    %remesh
    [x2,y2,tParametric,PARAM_2] = remesh(t,x2,y2,tParametric,PARAM_2);
    
    %current initial condition
    x0 = x2;
    y0 = y2;
    PARAM_0 = PARAM_2;

    %compute solution at t0
    [u0,v0,yLaplace,yStokes] = f(t,x0,y0,PARAM_0);
    
    if i~=loop
        dt = manyDT(i);

        %advance half step
        [x1,y1,PARAM_1] = advanceBlocks(x0,y0,tParametric,dt/2,u0,v0,PARAM_0);
        
        %compute solution at half step
        [u1,v1] = f(t+dt/2,x1,y1,PARAM_0);
        
%         hold off
%         plotGeometryDrop(x1,y1,PARAM_1,0)
%         drawnow

        %advance x0 with u1
        [x2,y2,PARAM_2] = advanceBlocks(x0,y0,tParametric,dt,u1,v1,PARAM_1);
        
    end
    
    %save data
    if sum((t<TTT+dt/10) & (t>TTT-dt/10))==1
        
        %staus function
        status = event(t,x0,y0,PARAM_0);
        if status==1
            %store results
            T(count) = t;
            X{count} = x0;
            U{count} = u0;
            Y{count} = y0;
            V{count} = v0;
            laplace{count} = yLaplace;
            stokes{count} = yStokes;
            T = T(1:count);
            X = X{1:count};
            U = U{1:count};
            Y = Y{1:count};
            V = V{1:count};
            break;
        end
        
        %store results
        T(count) = t;
        X{count} = x0;
        U{count} = u0;
        Y{count} = y0;
        V{count} = v0;
        laplace{count} = yLaplace;
        stokes{count} = yStokes;
        
        %increase counter
        count = count+1;
        
    end

end