

function [xNew,yNew] = splitElementsBubble(x,y,panel1,distance,PARAM)

count = 1;
xNew = x;
yNew = y;
dx = diff(x);   dy = diff(y);   dl = sqrt(dx.^2+dy.^2);
l = [0 cumsum(dl)];
ppx = spline(l,x);
ppy = spline(l,y);

removeNow = 0;
for i = 1:numel(x)-1
        
       %element coordinates
       x1 = x(i);
       x2 = x(i+1);
       y1 = y(i);
       y2 = y(i+1);
       l1 = l(i);
       l2 = l(i+1);
       lm = (l1+l2)/2;
       dl = sqrt((x1-x2)^2+(y1-y2)^2);

       %add node
       if (dl>distance(i)/PARAM.coeffDist) && (dl>PARAM.minSizeElemRemesh(panel1))
           
           display('Add new node')
           
           xNew(count) = x(i);
           yNew(count) = y(i);
           
           %build new node
           xNew(count+1) = ppval(ppx,lm);
           yNew(count+1) = ppval(ppy,lm);
        
            %new node position
            count = count+1;
        
       %remove node
%        elseif (dl<distance(i)/PARAM.coeffDist/2) && (dl<PARAM.maxElem(panel1)/2)
%            
%            display('Remove node')
%            
%            %don't remove
%            if removeNow==1
%                
%                xNew(count) = x(i);
%                yNew(count) = y(i);
%                removeNow = 0;
%                
%            else
%            
%                 count = count-1;
%                 %i = i+1;
%                 removeNow = 1;
%            
%            end
           
       else
           
           xNew(count) = x(i);
           yNew(count) = y(i);

       end
       count = count+1;
end

if xNew(end)~=x(end)
    xNew(numel(xNew)+1) = x(end);
    yNew(numel(yNew)+1) = y(end);
end
