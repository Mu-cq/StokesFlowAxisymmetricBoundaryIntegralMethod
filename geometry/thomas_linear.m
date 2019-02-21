function [ x ] = thomas_linear( a,b,c,s )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
K = length(a);
d = a;
y = a;
d(1) = 1/a(1)*b(1);
y(1) = 1/a(1)*s(1);

for i=1:(K-1)
    coef = 1./(a(i+1)-c(i)*d(i));
    if(i~=(K-1))
    d(i+1) = coef*b(i+1);
    end
    y(i+1) = coef*(s(i+1)-c(i)*y(i));
end

x(K) = y(K);
for i=(K-1):-1:1
    x(i) = y(i) - d(i)*x(i+1);
end
end

