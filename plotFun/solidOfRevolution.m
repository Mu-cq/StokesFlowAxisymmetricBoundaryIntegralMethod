%plot solid of revolution

function solidOfRevolution(r,z,theta,color,transp)

nGrid = 100;

%theta start end end
t1 = theta(1);
t2 = theta(2);

theta = linspace(t1,t2,nGrid);
theta = repmat(theta',1,numel(z));
z = repmat(z,nGrid,1);
r = repmat(r,nGrid,1);

%rotate shape
x = r.*cos(theta);
y = r.*sin(theta);

%plot3(x,y,z,'k')
%surf(z,y,x,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
%surf(z,y,x,'FaceColor',[1 1 1],'EdgeColor','none','LineStyle','none','FaceLighting','phong')
if color==10
    
    h = surf(z,y,x,zeros(size(x)),'FaceLighting','gouraud','LineStyle','none','AmbientStrength',0.5);
    
else
    
    h = surf(z,y,x,'FaceColor',color,'FaceLighting','gouraud','LineStyle','none','AmbientStrength',0.5);

end

%light('Position',[5 -5 5],'Style','inifinite')
%light('Position',[5 -5 5])
light('Position',[-10 -10 10])
material metal
alpha(h,transp)







