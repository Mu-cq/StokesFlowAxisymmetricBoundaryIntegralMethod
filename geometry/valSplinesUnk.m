%compute first and second derivative of a function computed with splines

function f = valSplinesUnk(t,a,b,c,d)

f = a+b.*t+c.*t.^2+d.*t.^3;