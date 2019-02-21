function stop = outfun(x,optimValues,state)
stop = false;
 
   switch state
       case 'init'
           hold on
       case 'iter'
          
           figure(4)
           if optimValues.iteration>0
           hold on
           end
           loglog(optimValues.iteration,optimValues.fval,'ok-')
           hold off
           grid on
           xlabel('ite')
           ylabel('fval')
           title('cost function value')
           drawnow
           
           figure(5)
           if optimValues.iteration>0
           hold on
           end
           loglog(optimValues.iteration,optimValues.firstorderopt,'or-')
           grid on
           xlabel('ite')
           ylabel('1st order opt')
           title('first order optimality condition')
           hold off
           drawnow
           
           %display(optimValues)
           
       case 'done'
           hold off
       otherwise
   end
end