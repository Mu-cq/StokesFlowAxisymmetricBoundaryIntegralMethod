%choose non linear function, 1 drop, Spectral

function [initial,x,y,newCa,xBase,yBase,perturb,PARAMupload] = initialConditionDrop(PARAM)

if PARAM.uploadShape==0
    
    if PARAM.algorithm==3||PARAM.algorithm==4
        
        disp('Build geometry from Acrivos 1973 linear theory')
        %build geometry from Acrivos 1973 linear theory
        Dguess = GetEllipseAcrivos(PARAM.Ca,PARAM.visc);
        
        if isnan(Dguess)
            error('Linear theory from Acrivos gives nan. Check input physical parameters')
        end
        
        %Dguess = 0.47;
        if PARAM.BC==2
            Dguess = PARAM.D;
        end
        [x,y] = ellipseCartesian(pi*PARAM.remeshMapping(PARAM.t),Dguess);
    
    elseif PARAM.algorithm==1||PARAM.algorithm==2||PARAM.algorithm==6

        disp('Draw initial shape analytically')

        if PARAM.shapeEllipse==1
            %build ellipse in curvilinear coordinates
            %semi axis from D and volume of sphere with radius=1
            D = PARAM.D;
            [x,y] = ellipseCartesian(pi*PARAM.remeshMapping(PARAM.t),D);
        elseif PARAM.shapeEllipse==0
            %build sphere perturb by legendre
            [x,y] = spherePlusLegendreCurvilinear(PARAM.remeshMapping(PARAM.t),PARAM.D,PARAM);
        elseif PARAM.shapeEllipse==2
            %build sphere perturb by 2 legendre polynimials P2 P4
            [x,y] = spherePlus2LegendreCurvilinear(PARAM.remeshMapping(PARAM.t),PARAM.D,PARAM.f2,PARAM);
        elseif PARAM.shapeEllipse==3
            %build sphere perturb by legendre P3, asymmetric
            [x,y] = spherePlusLegendreAsymmetric(PARAM.remeshMapping(PARAM.t),PARAM.D,PARAM);
        elseif PARAM.shapeEllipse==4
            %build sphere perturb by 2 legendre polynimials P2 P3
            [x,y] = spherePlus2LegendreCurvilinearAsymmetric(PARAM.remeshMapping(PARAM.t),PARAM.D,PARAM.f2,PARAM);
        end
        
    elseif PARAM.algorithm==5
    
        error('For stability analysis you have to upload a base state')
        
    else
        
        error('Missing option')
        
    end
    
    %dummy
    newCa = [];
    xBase = []; yBase = []; perturb = [];   PARAMupload = [];

elseif PARAM.uploadShape==1 || PARAM.uploadShape==2
    
    if PARAM.uploadShape==1
    
        disp('Upload initial shape from continuation')

        %upload shape
        [x,y,newCa,xBase,yBase,perturb,PARAMupload] = uploadFromContinuation(PARAM);
        
    elseif PARAM.uploadShape==2
        
        disp('Upload initial shape from single newton method')
        
        %upload shape
        [x,y,xBase,yBase,perturb,PARAMupload] = uploadFromNewton(PARAM);
        
        %dummy
        newCa = PARAM.Ca;
        
    end
    
    if PARAM.overlapMode==1
            %perturb Shape with Legendre polynomial
            [x,y] = perturbWithLegendrePoly(x,y,PARAM.D,PARAM);
    elseif PARAM.overlapMode==2
            %perturb Shape with eigenmode
            if PARAM.uploadShape==1
                error('Overlap-mode works well only when uploading from one Newton Method')
            end
            [x,y] = perturbEigenmode(x,y,PARAM.D,PARAM);
    end
    
elseif PARAM.uploadShape==3
    
    disp('Upload shape from edge tracking')
    
    %upload shape
    [x,y,xBase,yBase,perturb,PARAMupload] = uploadFromEdgeTracking(PARAM);
    
    if PARAM.overlapMode==1
            %perturb Shape with Legendre polynomial
            [x,y] = perturbWithLegendrePoly(x,y,PARAM.D,PARAM);
    elseif PARAM.overlapMode==2
            %perturb Shape with eigenmode
            error('This works well only when uploading from one Newton Method')
    end
    
    %dummy
    newCa = PARAM.Ca;
    
elseif PARAM.uploadShape==4
    
    disp('Upload shape from DNS')
    
    %upload shape
    [x,y,xBase,yBase,perturb,PARAMupload] = uploadFromDNS(PARAM);
    
    if PARAM.overlapMode==1
            %perturb Shape with Legendre polynomial
            [x,y] = perturbWithLegendrePoly(x,y,PARAM.D,PARAM);
    elseif PARAM.overlapMode==2
            %perturb Shape with eigenmode
            error('This works well only when uploading from one Newton Method')
    end
    
    %dummy
    newCa = PARAM.Ca;
    
elseif PARAM.uploadShape==5
    
    disp(['Upload shape from Minimal Seed: A0=' num2str(PARAM.A0upload) ' and T=' num2str(PARAM.ThorizonUpload)])
    
    %upload shape
    [x,y,xBase,yBase,perturb,PARAMupload] = uploadFromMinimalSeed(PARAM);
    
end

%place shape in the center of mass
if PARAM.placeShapeXCM==1
    
    disp('Place shape center of mass in the origin')
    xcm = CenterMassCurvAxisSpectral(x,y,PARAM);
    x = x-xcm;
    
end

%go from grid points to modes
[xMode,yMode] = fromGridToModes(x,y,PARAM);
xMode(PARAM.dealiasing+1:end) = 0;  yMode(PARAM.dealiasing+1:end) = 0;
initial = zeros(2*PARAM.n+2,1);
initial(1:2:end-1) = xMode;     initial(2:2:end) = yMode;

if PARAM.remeshStart==1
    %execute remesh
    disp('Remesh spectrally')
    initial = remeshDropSpectral(initial,PARAM);
end

