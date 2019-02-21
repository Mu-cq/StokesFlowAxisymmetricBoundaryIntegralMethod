%clear all
%close all

function rOut = findRbubble(Ad2,rhs,in,rBubble)

eqn = Ad2*rBubble^3 + 2*rBubble^2 - rhs == 0;

sol = vpasolve(eqn,in);

if imag(eval(sol(1,1)))==0 && real(eval(sol(1,1)))>0
    rOut = eval(sol(1,1));
elseif imag(eval(sol(2,1)))==0 && real(eval(sol(2,1)))>0
    rOut = eval(sol(2,1));
elseif imag(eval(sol(3,1)))==0 && real(eval(sol(3,1)))>0
    rOut = eval(sol(3,1));
else
    rOut = 0;
end

end