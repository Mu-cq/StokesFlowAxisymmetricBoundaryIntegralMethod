%draw bubble in channel

function [x,y] = drawBubbleInChannel(x0,y0,rBubble,hFilm,Hchannel,nBubble)

if rBubble<Hchannel-hFilm
    
   [x,y] = draw_circle_clean(x0,y0,nBubble,0);
   x = x*rBubble;
   y = y*rBubble;
    
else
    
   %draw squeezed bubble
   V = 4/3*pi*rBubble^3;
   R = Hchannel-hFilm;
   L = (V-4/3*pi*R^2)/pi/R^2;
    
   theta = linspace(0,pi/2,100);
   L2 = linspace(0,1,100)*L;
   L1 = R*theta;    L3 = L1;
   Lcurv = [L1 L2(2:end)+L1(end) L3(2:end)+L2(end)+L1(end)];
   
   xTemp = [R*cos(theta)+L/2+x0 L/2-L2(2:end)+x0 R*cos(theta(2:end)+pi/2)-L/2+x0];
   yTemp = [R*sin(theta)+y0 (Hchannel-hFilm)*ones(1,numel(L2(2:end)))+y0 R*sin(theta(2:end)+pi/2)+y0];
   
   x = spline(Lcurv,xTemp,linspace(0,Lcurv(end),nBubble+1));
   y = spline(Lcurv,yTemp,linspace(0,Lcurv(end),nBubble+1));
   
end