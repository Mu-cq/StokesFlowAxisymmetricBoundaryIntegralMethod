%clear all
%close all

function [R1,R2,H] = findR1R2Hmod(PARAM)

%solve system of algebraic equation

%input data
theta = PARAM.theta;
h = PARAM.h;
L = PARAM.L;
R = PARAM.R;
Vtot = PARAM.Vtot;
alpha1 = pi/2+theta;
alpha2 = pi/2-theta;
x0 = PARAM.start;

%syms R1 R2 H

eqn1 = (R1+h)*cos(theta)-tan(theta)*(x0+H/2-(R1+h)*sin(theta)-L/2) - R == 0;
eqn2 = (R2+h)*cos(theta)-tan(theta)*(x0-H/2-(R2+h)*sin(theta)-L/2) - R == 0;
eqn3 = 1/3*pi*(H+R2*sin(theta)-R1*sin(theta))*(R1^2*cos(theta)^2+R2^2*cos(theta)^2+R1*R2*cos(theta)^2)...
    - pi/3*cos(theta)^2*sin(theta)*(R2^3-R1^3)...
    +2*pi/3*R1^3*(-cos(alpha1)+1) +2*pi/3*R2^3*(-cos(alpha2)+1) == Vtot;

sol = vpasolve(eqn1, eqn2, eqn3, [2 1 1]);

tempH = eval(sol.H);
tempR1 = eval(sol.R1);
tempR2 = eval(sol.R2);

H = real(tempH(find(abs(imag(tempH))<10e-8==1),:));
R1 = real(tempR1(find(abs(imag(tempR1))<10e-8==1),:));
R2 = real(tempR2(find(abs(imag(tempR2))<10e-8==1),:));

end