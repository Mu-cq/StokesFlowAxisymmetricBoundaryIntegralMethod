%get shape based on some norm criteria and averaged them applying volume
%conservation

function [Tstart,initialShape] = getShapesAndAverageBEM(T1,Y1,T2,Y2,PARAM,V0)

Tcell = cell(2,1);
lNorm = cell(2,1);
Tcell{1} = T1;  Tcell{2} = T2;
Ycell{1} = Y1;  Ycell{2} = Y2;

%compute norm of the shapes
for k = 1:2
    
     T = Tcell{k};
     Y = Ycell{k};
     L = zeros(numel(T),1);
     for i = 1:numel(T)
            
         %load shape
         Yhere = Y{i};
         x = Yhere(1:2:end-1)';
         y = Yhere(2:2:end)';
            
         %radius minus sphere
         xcm = center_mass(x,y);
         f = sqrt((x-xcm).^2+y.^2)-1;
         V = axis_int_gauss_vect(x,y);
         R0 = nthroot(V/4/pi*3,3);
         if PARAM.normEdge==1    %L-2 norm
                weight = integrationOnLineWeightAxis(x,y,PARAM.orderVariableStokes(PARAM.blockEdge),PARAM.orderGeometryStokes(PARAM.blockEdge),PARAM.SPlinesType(PARAM.blockEdge));
                L(i) = sqrt(weight*(f.^2));
         elseif PARAM.normEdge==2    %L-infinity norm of radius relative to center of mass
                L(i) = max(abs(f));
         elseif PARAM.normEdge==3    %excess surface norm
                L(i) = surf_gauss_vect(x,y)-4*pi*R0^2;
         end
        
     end
     lNorm{k} = L;
    
end

lnorm1 = lNorm{1};
lnorm2 = lNorm{2};

BreakNow = 0;
firstShape = 0;
for i = 1:numel(lnorm1)
    
    for k = 1:numel(lnorm2)
        
        if abs(T1(i)-T2(k))<1e-6
            
            %get the default value for bi-section
            if firstShape==0
                
                NewShape1 = i;
                NewShape2 = k;
                
            end
            firstShape = 1;
            
            %difference between the norms
            normDiff = abs(lnorm1(i)-lnorm2(k));

            %check if norm is larger
            if normDiff>PARAM.deltaEdge
                if isempty(PARAM.forceIntegrationAfterDT)==1 || T1(i)==T1(1) || T2(i)==T2(1)
                    disp(['New shape is taken at T=' num2str(T1(i)) ', where Norm2-Norm1=' num2str(normDiff)])
                end
                NewShape1 = i;
                NewShape2 = k;
                BreakNow = 1;
                break;
            end
        
        end
    
    end
    
    if BreakNow==1
        break;
    end
    
end

%force to bisection after a fixed time DT
if isempty(PARAM.forceIntegrationAfterDT)==0 && NewShape1>1 && NewShape2>1
    
    newTime = max(T1(1),T2(1))+PARAM.forceIntegrationAfterDT;
    
    disp(['The bisection tolerance has been reached, next shape is taken at T=' num2str(newTime)])
    
    [~,NewShape1] = min(abs(newTime-T1));
    [~,NewShape2] = min(abs(newTime-T2));
    
end

%select shape
Y1here = Y1{NewShape1};
Y2here = Y2{NewShape2};
x1 = Y1here(1:2:end-1)';
y1 = Y1here(2:2:end)';
x2 = Y2here(1:2:end-1)';
y2 = Y2here(2:2:end)';

%place shape in the center of mass
if PARAM.placeShapeXCM==1
   
   disp('Place shape (to bisect) center of mass in the origin')
   xcm1 = center_mass(x1,y1);
   xcm2 = center_mass(x2,y2);
   x1 = x1-xcm1;
   x2 = x2-xcm2;
    
end

%average shapes
if numel(x1)~=numel(x2)
    
    [x1,y1,x2,y2] = equalNumberOfNodes(x1,y1,x2,y2,1,PARAM);
    
end
x = (x1+x2)/2;
y = (y1+y2)/2;

%compute normal vector
[nx,ny] = computeNormalVector(x,y,PARAM.orderVariableStokes(PARAM.blockEdge),PARAM.orderGeometryStokes(PARAM.blockEdge),PARAM.SPlinesType(PARAM.blockEdge));

if PARAM.bisection==1
    
    disp('Bisect by simple average')

    %compute modes first y mode from volume (x is set to zero)
    fVolume = @(rho) ModifyVolumeDNS(x,y,nx,ny,rho,V0);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    displ = fsolve(fVolume,0,options);

    %displace in the normal direction
    x = x + displ*nx;
    y = y + displ*ny;
    
elseif PARAM.bisection==2
    
    disp('Bisect area')
    
    %desired area
    A1 = surf_gauss_vect(x1,y1);
    A2 = surf_gauss_vect(x2,y2);
    Afinal = (A1+A2)/2;
    
    %second Legendre mode
    theta = atan(y./x);
    theta = theta + pi*(theta<0);
    P2 = legendre(2,cos(theta));
    P2 = P2(1,:);
    
    %compute rho in symmetry axis
    fVolume = @(unk) ModifyVolumeAreaXY(x,y,nx,ny,unk,V0,Afinal,P2);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    move = fsolve(fVolume,[0 0],options);

    %compute full shape
    x = x + nx*move(1) + nx*move(2).*P2;
    y = y + ny*move(1) + ny*move(2).*P2 ;
    
end

%output
Tstart = T1(NewShape1);
initialShape = zeros(2*numel(x),1);
initialShape(1:2:end-1) = x;
initialShape(2:2:end) = y;








