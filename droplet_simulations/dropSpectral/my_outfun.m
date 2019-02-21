% my output function for fsolve

function stop = my_outfun(x,optimValues,state)

stop = 0;

% compute residuals
res = norm(optimValues.fval, Inf);

% plot residuals
if strcmp(state,'init') || strcmp(state,'iter')
    semilogy(optimValues.iteration,res,'ok')
    hold on
elseif strcmp(state,'done')
    grid on
    hold on
    xyLabelTex('\rm{ite}','||\hat{\mathbf{u}}^{(n)}||_\infty')
    title('Residuals')
end