%group solution by panel

function xNew = groupSolutionByPanel(x,nnn)

xNew = cell(numel(nnn),1);
for i = 1:numel(nnn)
    
    if i==1
          startMatrix = 1;
    else
          startMatrix = 1+sum(nnn(1:i-1))*2;
    end
    endMatrix = sum(nnn(1:i))*2;
    
    xNew{i} = x(startMatrix:endMatrix)';
    
end