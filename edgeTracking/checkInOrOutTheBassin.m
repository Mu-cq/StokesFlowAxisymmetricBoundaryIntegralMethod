%get shape based on some norm criteria and averaged them applying volume
%conservation

function [T1out,Y1out,T2out,Y2out] = checkInOrOutTheBassin(Tcell,Ycell,PARAM)

DDD = zeros(numel(Tcell),1);

%compute norm of the shapes
for k = 1:numel(Tcell)
    
     T = Tcell{k};
     Y = Ycell{k};
     i = numel(T);
            
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
            
     %defromation parameter at last shape
     L = max(x)-min(x);
     B = 2*max(y);
     D = (L-B)/(L+B);
     DDD(k) = D;
    
end

%check if I'm inside or ouside the bassin of attraction
if PARAM.BC==1
    
    if PARAM.visc==1
        load('./steadyState/CaExt')
        load('./steadyState/DExt')
    elseif PARAM.visc==0
        load('./steadyState/CaExt0')
        load('./steadyState/DExt0')
    elseif PARAM.visc==0.1
        load('./steadyState/CaExt01')
        load('./steadyState/DExt01')
    elseif PARAM.visc==0.5
        load('./steadyState/CaExt05')
        load('./steadyState/DExt05')
    elseif PARAM.visc==5
        load('./steadyState/CaExt5')
        load('./steadyState/DExt5')
    elseif PARAM.visc==10
        load('./steadyState/CaExt10')
        load('./steadyState/DExt10')
    end
    [~,ind] = min(abs(PARAM.Ca-manyCa));
    Dca = manyD(ind);
    OutIn  = abs(DDD-Dca)/Dca > PARAM.convergeShapeEdge;
    
elseif PARAM.BC==2
    
    OutIn  = abs(DDD) > PARAM.convergeShapeEdge;
   
end

if sum(OutIn==1)==numel(OutIn)
    error('All shapes are outside the bassin of attraction')
elseif sum(OutIn==1)==0
    error('All shapes are inside the bassin of attraction')
end

%fin last stable
indIN = find(OutIn==0,1,'last');

%fin last unstable
indOUT = find(OutIn==1,1,'last');

display(['Average between stable shape ' num2str(indIN) ' and unstable shape ' num2str(indOUT)])

%output
T1out = Tcell{indIN};
Y1out = Ycell{indIN};
T2out = Tcell{indOUT};
Y2out = Ycell{indOUT};







