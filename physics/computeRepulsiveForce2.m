%compute Repulsive Force

function [F,dist] = computeRepulsiveForce2(x,y,block1,block2,PARAM)

%error('Still to be debugged')

%get panels for compting the distance
[xBlock1,yBlock1,~,~,manyPanel1] = getBlockCoordinates(x,y,PARAM,block1);
[xBlock2,yBlock2,~,~,manyPanel2] = getBlockCoordinates(x,y,PARAM,block2);

%find closet panels
minTemp = zeros(numel(manyPanel2),1);
minHere = zeros(numel(manyPanel1),1);
indHere = zeros(numel(manyPanel1),1);
for i = 1:numel(manyPanel1)
    
    for k = 1:numel(manyPanel2)
        
        %minTemp(k) = panelDistance(x,y,manyPanel1(i),manyPanel2(k),PARAM);
        x1 = x{manyPanel1(i)};  y1 = y{manyPanel1(i)};
        x2 = x{manyPanel2(k)};  y2 = y{manyPanel2(k)};
        x1 = (x1(1:end-1)+x1(2:end))/2; y1 = (y1(1:end-1)+y1(2:end))/2;
        x2 = (x2(1:end-1)+x2(2:end))/2; y2 = (y2(1:end-1)+y2(2:end))/2;
        minTemp(k) = min(distWallDrop(x1,y1,x2,y2));
        
    end
    
    [minHere(i),indHere(i)] = min(minTemp);
    
end
[~,ind] = min(minHere);
panel1 = manyPanel1(ind);
panel2 = manyPanel2(indHere(ind));

%compute smallest separation
[dist,~,~,sign] = panelDistance(x,y,panel1,panel2,PARAM);
%[dist,sign] = geometryRepulsive(xBlock1,yBlock1,xBlock2,yBlock2);%   sign = -sign;

%critical distance
delta = PARAM.repulsiveOn;
r = dist;

%compute only when surfaces are close
activate = r<=delta;

%compute intensity repulsive force
F = PARAM.coeffRepulsive*(exp(delta-r)-1)/(exp(delta)-1)*activate*sign(1)/abs(sign(1));



