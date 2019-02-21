%analytical solution of poisson eqution in 2D

clear all
close all

%domain
L = 10;
R = 1;
x = linspace(0,L,10*L);
y = linspace(0,R,10*R);
[X,Y] = meshgrid(x,y);
a = 5; b = 0; %place the dirac

%numeber of modes
nk = 101;
nj = 101;

phi = 0;
%compution with fourier series
for i = nk:-1:1
    
    %Ak = sin((i-1)*pi*b/R);
        
    for l = nj:-1:1
        
        %Bj = cos((l-1)*pi*b/L);
                
        %Ckj = Ak*Bj;
        Ckj = sin((l-1)*pi*a/L)*cos((i-1)*pi*b/R);
        phi = phi + Ckj*cos(pi*(i-1)*Y/R).*sin(pi*(l-1)*X/L);
    
    end
end

figure
contourf(X,Y,phi)
axis equal