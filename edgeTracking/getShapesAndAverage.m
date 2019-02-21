%get shape based on some norm criteria and averaged them applying volume
%conservation

function [Tstart,initialShape] = getShapesAndAverage(T1,Y1,T2,Y2,PARAM,V0)

Tcell = cell(2,1);
lNorm = cell(2,1);
Tcell{1} = T1;  Tcell{2} = T2;
Ycell{1} = Y1;  Ycell{2} = Y2;
w = PARAM.WG';

%compute norm of the shapes
for k = 1:2
    
     T = Tcell{k};
     Y = Ycell{k};
     L = zeros(numel(T),1);
     for i = 1:numel(T)
            
         %load modes
         xMode = Y(i,1:2:end-1)';
         yMode = Y(i,2:2:end)';

         %build shape
         if PARAM.legendre==1

                   x = LegendreBuildXY(xMode,PARAM.PPP);
                   y = LegendreBuildXY(yMode,PARAM.PPP);

         elseif PARAM.legendre==0

                   x = chebcoeffs2chebvals(xMode);
                   y = chebcoeffs2chebvals(yMode);

         end
            
         %radius minus sphere
         xcm = CenterMassCurvAxisSpectral(x,y,PARAM);
         f = sqrt((x-xcm).^2+y.^2)-1;
         V = VolumeCurvilinearAxisSpectral(x,y,PARAM);
         R0 = nthroot(V/4/pi*3,3);
         if PARAM.normEdge==1    %L-2 norm
                L(i) = sqrt(w*(f.^2));
         elseif PARAM.normEdge==2    %L-infinity norm
                L(i) = max(abs(f));
         elseif PARAM.normEdge==3    %excess surface norm
                L(i) = surfaceCurvilinearAxisSpectral(x,y,PARAM)-4*pi*R0^2;
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
            
            if PARAM.bruteForceBisection==1
                display(['The bisection criteria has been by-passed, Norm2-Norm1=' num2str(normDiff)])
                NewShape1 = i;
                NewShape2 = k;
                BreakNow = 1;
                break;
            end

            %check if norm is larger
            if normDiff>PARAM.deltaEdge*PARAM.deltaModify
                fprintf('New shape is taken at T=%f, where Norm2-Norm1=%1.16f \n',T1(i),normDiff)
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

%select shape
xMode1 = Y1(NewShape1,1:2:end-1)';
yMode1 = Y1(NewShape1,2:2:end)';
xMode2 = Y2(NewShape2,1:2:end-1)';
yMode2 = Y2(NewShape2,2:2:end)';

%place shape in the center of mass
if PARAM.placeShapeXCM==1
   
   disp('Place shape center of mass in the origin')
   [x1,y1] = fromModesToGrid(xMode1,yMode1,PARAM);
   [x2,y2] = fromModesToGrid(xMode2,yMode2,PARAM);
   xcm1 = CenterMassCurvAxisSpectral(x1,y1,PARAM);
   xcm2 = CenterMassCurvAxisSpectral(x2,y2,PARAM);
   x1 = x1-xcm1;
   x2 = x2-xcm2;
   [xMode1,yMode1] = fromGridToModes(x1,y1,PARAM);
   [xMode2,yMode2] = fromGridToModes(x2,y2,PARAM);
    
end

%average modes
xMode = (xMode1+xMode2)/2;
yMode = (yMode1+yMode2)/2;

%compute grid points
[x,y] = fromModesToGrid(xMode,yMode,PARAM);

%compute normal vector
[nx,ny] = normalVectorSpectral(x,y,PARAM);

if PARAM.bisection==1
    
    disp('Bisect by simple average')

    %compute modes first y mode from volume (x is set to zero)
    fVolume = @(rho) ModifyVolumeSpectralXYdns(x,y,nx,ny,rho,V0,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    displ = fsolve(fVolume,0,options);

    %displace in the normal direction
    x = x + displ*nx;
    y = y + displ*ny;
    
elseif PARAM.bisection==2
    
    disp('Bisect area')
    
    %modes
    if PARAM.legendre==1||PARAM.legendre==2
        PPP = PARAM.PPP;
    elseif PARAM.legendre==0
        PPP = PARAM.TTT;
    end
    
    %desired area
    [x1,y1] = fromModesToGrid(xMode1,yMode1,PARAM);
    [x2,y2] = fromModesToGrid(xMode2,yMode2,PARAM);
    A1 = surfaceCurvilinearAxisSpectral(x1,y1,PARAM);
    A2 = surfaceCurvilinearAxisSpectral(x2,y2,PARAM);
    Afinal = (A1+A2)/2;
    
    %compute rho in symmetry axis
    fVolume = @(unk) ModifyVolumeAreaSpectralXY(x,y,nx,ny,unk,V0,Afinal,PARAM);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    move = fsolve(fVolume,[0 0],options);

    %compute full shape
    x = x + nx*move(1).*PPP(1,:)' + nx*move(2).*PPP(3,:)';
    y = y + ny*move(1).*PPP(1,:)' + ny*move(2).*PPP(3,:)';
    
end

%place first and last node on the axis
% if PARAM.legendre==0 || PARAM.legendre==2
%     disp('Impose that first and last node are on the axis')
%     y([1 end]) = [0 0];
% end
    
%dealiasing
[x,y] = dealiasingGridXY(x,y,PARAM);

%from gird to points modes
[xMode,yMode] = fromGridToModes(x,y,PARAM);

%output
Tstart = T1(NewShape1);
initialShape = zeros(2*numel(xMode),1);
initialShape(1:2:end-1) = xMode;
initialShape(2:2:end) = yMode;








