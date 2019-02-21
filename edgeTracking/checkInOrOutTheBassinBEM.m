%get shape based on some norm criteria and averaged them applying volume
%conservation

function [T1out,Y1out,T2out,Y2out] = checkInOrOutTheBassinBEM(Tcell,Ycell,PARAM)

DDD = zeros(numel(Tcell),1);

%compute norm of the shapes
for k = 1:numel(Tcell)
    
     T = Tcell{k};
     Y = Ycell{k};
     i = numel(T);
            
     %load shape
     Y = Y{i};
     x = Y(1:2:end-1)';
     y = Y(2:2:end)';
            
     %defromation parameter at last shape
     L = max(x)-min(x);
     B = 2*max(y);
     D = (L-B)/(L+B);
     DDD(k) = D;
    
end

%check if I'm inside or ouside the bassin of attraction
if PARAM.typeBCstokes(PARAM.blockEdge)==2
    
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
    
elseif PARAM.typeBCstokes(PARAM.blockEdge)==7
    
    OutIn  = abs(DDD) > PARAM.convergeShapeEdge;
    
else
    
    error('Not implemented')
   
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

if PARAM.interpOrExtrap==1

    display(['Average between stable shape ' num2str(indIN) ' and unstable shape ' num2str(indOUT)])

elseif PARAM.interpOrExtrap==2
    
    display(['Extrapolate from shape ' num2str(indIN) ' and bisect one scalar parameter with unstable shape ' num2str(indOUT)])

else

    error('Not implemented')
    
end

%output
T1out = Tcell{indIN};
Y1out = Ycell{indIN};
T2out = Tcell{indOUT};
Y2out = Ycell{indOUT};







