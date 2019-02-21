%euqtion of motion for a bubble in a cone fue to lubrication equation

clear variables
close all

%parameter
theta = 0.1;

%motor geometry
H = 1;
L = 10;

%initial condition
xb = 11;
xm = 0;
a = 1;
ap = 1;

%time stepping
dt = 1e-3;
loop = 400;

%initialize
manyUB = zeros(loop,1);
manyUM = zeros(loop,1);
manyXB = zeros(loop+1,1);
manyXM = zeros(loop+1,1);
manyA = zeros(loop+1,1);

%initial condition
manyXB(1) = xb;
manyXM(1) = xm;
manyA(1) = a;

%euler loop
for i = 1:loop
    
    %display(['Loop ' num2str(i) ' of ' num2str(loop)])
    
    %compute gap
    gap = (xb-xm)*sin(theta)-a;

    if gap<0
        error('Gap size must be positive')
    end
    
    %compute dirty bubble and motor velocity
    [Ub,Um] = dirtyBubbleAndMotorVelocityLubrication(a,ap,theta,gap,L,H);
    
    %advance in time
    xb = xb + Ub*dt;
    a = a + ap*dt;
    xm = xm + Um*dt;
    
    %save data
    manyUB(i) = Ub;
    manyUM(i) = Um;
    manyXB(i+1) = xb;
    manyXM(i+1) = xm;
    manyA(i+1) = a;
    
end

%plot initial and final condition
figure
subplot(2,2,1)
thetaSphere = linspace(0,2*pi,100);
xCone = linspace(manyXB(1)-2,manyXB(end)+2,100);
yCone = xCone*tan(theta);
plot(manyA(1)*cos(thetaSphere)+manyXB(1),manyA(1)*sin(thetaSphere))
hold on
plot(manyA(end)*cos(thetaSphere)+manyXB(end),manyA(end)*sin(thetaSphere))
plot(xCone,yCone,'k')
plot(xCone,-yCone,'k')
plot(xCone,yCone-(xm*tan(theta)),'r')
%legend('Initial Shape','Final Shape','Location','Best')
axis equal
grid on
xlabel('z')
ylabel('r')
title(['\theta=' num2str(theta) ' z_b^{in}=' num2str(manyXB(1)) ' a^{in}=' num2str(manyA(1))])

%plot bubble velocity
%figure
subplot(2,2,2)
plot(0:dt:(loop-1)*dt,manyUB)
grid on
xlabel('time')
ylabel('U_b')
title('Bubble velocity')

%plot bubble velocity
%figure
subplot(2,2,3)
plot(0:dt:(loop-1)*dt,manyUM)
grid on
xlabel('time')
ylabel('U_m')
title('Motor velocity')

% %plot bubble position
% %figure
% subplot(2,2,3)
% plot(0:dt:loop*dt,manyXB)
% grid on
% xlabel('time')
% ylabel('z_b')
% title('Bubble position')

%plot gap size
%figure
subplot(2,2,4)
semilogy(0:dt:loop*dt,(manyXB-manyXM)*sin(theta)-manyA)
grid on
xlabel('time')
ylabel('\epsilon')
title('Gap size')









