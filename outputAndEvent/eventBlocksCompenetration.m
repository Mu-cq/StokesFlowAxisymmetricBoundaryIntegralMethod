%activate spectral remesh

function [value,isterminal,direction] = eventBlocksCompenetration(t,x,y,PARAM)

%check compenetration
PARAM.typeBC = PARAM.typeBCstokes;
PARAM.orderVariable = PARAM.orderVariableStokes;
PARAM.orderGeometry = PARAM.orderGeometryStokes;
value = 0;
for i = 1:numel(PARAM.blockType)
    
    [xHere1,yHere1] = getBlockCoordinates(x,y,PARAM,i);
    
    for k = 1:numel(PARAM.blockType)

        if k~=i
            [xHere2,yHere2] = getBlockCoordinates(x,y,PARAM,k);
            InOut = FigInOut(xHere1,yHere1,xHere2',yHere2');
            if max(InOut)>pi
                warning('Compenetration occurs')
                value = value+1;
            end
        end

    end

end

%count elements
elem = 0;
for i = 1:numel(PARAM.n)
    
    elem = elem+numel(x{i});
    
end

if value==0
    
    %output to screen
    fprintf('T=%f elem=%i\n',t,elem);

end

value = (value>0);
isterminal = 1;
direction = 0;