%advance blocks position

function [x,y,PARAM] = advanceBlocks(x,y,tParametric,dt,u,v,PARAM)

%loop on blocks
for i = 1:numel(PARAM.n)
    
   if PARAM.panelType(i)==1 %if it is a moving block
       
       if PARAM.panelInflate(i)==0
           
           if PARAM.geometryPanel(i)==0 %straight line
               
               PARAM.xStart(i) = PARAM.xStart(i) + dt*u{i};
               PARAM.xEnd(i) = PARAM.xEnd(i) + dt*u{i};
               PARAM.yStart(i) = PARAM.yStart(i) + dt*v{i};
               PARAM.yEnd(i) = PARAM.yEnd(i) + dt*v{i};
               
               x{i} = x{i}+dt*u{i};
               y{i} = y{i}+dt*v{i};
               
           elseif PARAM.geometryPanel(i)==1 %circle
                   
               PARAM.x0_Circle(i) = PARAM.x0_Circle(i) + dt*u{i};
               PARAM.y0_Circle(i) = PARAM.y0_Circle(i) + dt*v{i};
               
               x{i} = x{i}+dt*u{i};
               y{i} = y{i}+dt*v{i};
                   
           end
           
       elseif PARAM.panelInflate(i)==1
           
           if PARAM.geometryPanel(i)==0 %straight line
               
               error('Not Implemented')
               
               PARAM.xStart(i) = PARAM.xStart(i) + dt*u{i};
               PARAM.xEnd(i) = PARAM.xEnd(i) + dt*u{i};
               PARAM.yStart(i) = PARAM.yStart(i) + dt*v{i};
               PARAM.yEnd(i) = PARAM.yEnd(i) + dt*v{i};
               
           elseif PARAM.geometryPanel(i)==1 %circle
               
               uHere = u{i};
                   
               PARAM.x0_Circle(i) = PARAM.x0_Circle(i) + dt*uHere(1);
               PARAM.y0_Circle(i) = PARAM.y0_Circle(i) + dt*v{i};
               
               PARAM.rArc(i) = PARAM.rArc(i) + dt*uHere(2);
               
               x{i} = x{i}+dt*uHere(1) + dt*uHere(2)*cos(tParametric{i});
               y{i} = y{i}+dt*v{i} + dt*uHere(2)*sin(tParametric{i});
                   
           end
           
       end
       
   elseif PARAM.panelType(i)==2 % if it is a droplet
       
       x{i} = x{i}+dt*u{i};
       y{i} = y{i}+dt*v{i};
       
   end
    
end