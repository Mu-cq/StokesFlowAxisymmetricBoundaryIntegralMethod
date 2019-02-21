close all;
clc; clear;

% Force Dipole parallel or perpendicular to wall
% Streamlines visualization 

prompt = 'Choose x component of force 1 = ';
F1(1) = input(prompt);

F1(2) = 0;

prompt = 'Choose z component of force 1 = ';
F1(3) = input(prompt);

F2 = - F1;

prompt = 'Choose distance between F1 and F2 = ';
d = input(prompt);

phi = atan2(F1(3),F1(1));

% choose size of box
lx = 20;
lz = 20;

% Choose where the center of the dipole is placed
c = [5,0,2];

x1(1) = (d/2)*cos(phi) + c(1);
x1(2) = 0;
x1(3) = (d/2)*sin(phi) + c(3);

x2(1) = -(d/2)*cos(phi) + c(1);
x2(2) = 0;
x2(3) = -(d/2)*sin(phi) + c(3);

x = 0:0.25:lx; z = 0:0.25:lz;
[X,Z] = meshgrid(x,z);

u = zeros(length(x),length(z));
v = zeros(length(x),length(z));
w = zeros(length(x),length(z));

deltatensor = eye(3,3);
deltatensor(3,3) = -1;

identity = eye(3,3);
deltai3 = [0,0,1];

% Choose viscosity
mu = 1;
% F1 = F1/(mu*pi*8);
% F2 = F2/(mu*pi*8);

F1delta = F1*deltatensor;
F2delta = F2*deltatensor;

for i = 1:1:length(x)
        for j = 1:1:length(z)
            
% compute the  flow field for every grid point for Stokeslet 1
            h = abs(x1(3));
            
            r(1) = X(i,j) - x1(1);
            r(2) = 0;
            r(3) = Z(i,j) - x1(3);
            
% Wall is located at z=0 
            
            R(1) = X(i,j) - x1(1);
            R(2) = 0;
            R(3) = Z(i,j) + x1(3);
                             
            rr = norm(r)+0.001;
            RR = norm(R);
                        
            U1 = F1/rr + dot(F1,r)*r/rr^3 - (F1/RR + dot(F1,R)*R/RR^3) ...
                + 2*h*h*(F1delta/RR^3 - 3*dot(F1delta,R)*R/RR^5) ...
                + 2*h*(deltai3*dot(F1delta,R)/RR^3 - F1delta*R(3)/RR^3 ...
                - F1delta(3)*R/RR^3 + 3*R*R(3)*dot(F1delta,R)/RR^5);
            
%             if rr<=0.1
%                 U1=[0,0,0];
%             end

            % Net velocity is 
            u(i,j) = U1(1);
            w(i,j) = U1(3);
% compute the  flow field for every grid point for Stokeslet 2

            h = abs(x2(3));

            r(1) = X(i,j) - x2(1);
            r(2) = 0;
            r(3) = Z(i,j) - x2(3);
%             
% Wall is located at z=0 

            R(1) = X(i,j) - x2(1);
            R(2) = 0;
            R(3) = Z(i,j) + x2(3);
                             
            rr = norm(r)+0.001;
            RR = norm(R);
            
            U2 = F2/rr + dot(F2,r)*r/rr^3 - (F2/RR + dot(F2,R)*R/RR^3) ...
                + 2*h*h*(F2delta/RR^3 - 3*dot(F2delta,R)*R/RR^5) ...
                + 2*h*(deltai3*dot(F2delta,R)/RR^3 - F2delta*R(3)/RR^3 ...
                - F2delta(3)*R/RR^3 + 3*R*R(3)*dot(F2delta,R)/RR^5);
%             
% %             if rr<=1
% %                 U2=[0,0,0];
% %             end
% 
%             % Net velocity is 
%             u(i,j) = U1(1) + U2(1);
%             w(i,j) = U1(3) + U2(3);


        end
        
        
end

h1 = figure(1); 


startx = x;
startz = zeros(size(startx));
hlines = streamline(x,z,u,w,startz,startx);
set(hlines,'LineWidth',1.5,'Color','blue');
   
% startz = zeros(size(x));
% hlines = streamline(x,z,u,w,x,startz);
% set(hlines,'LineWidth',1.5,'Color','blue');
    
% startz = 0*ones(size(x));
% hlines = streamline(x,z,u,w,x,startz);
% set(hlines,'LineWidth',1.5,'Color','blue')
% 
startz = 20*ones(size(startx));
hlines = streamline(x,z,u,w,startz,startx);
set(hlines,'LineWidth',1.5,'Color','blue')
% 
% startz = 10*ones(size(x));
% hlines = streamline(x,z,u,w,x,startz);
% set(hlines,'LineWidth',1.5,'Color','blue')

%streamslice(x,z,u,w)

axis equal;
axis([0 10 0 10 0 10]);
view([0 0 1])
set(gca,'xtick',[],'ytick',[]);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

figuresize(18,18,'centimeters');
saveas(h1,'forcedipole3','pdf');




















