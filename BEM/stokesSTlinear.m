%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [SXX,SYY,anA,anB] = stokesSTlinear(SXX,SYY,dL,x,y,x0,y0)

    % number of singularities on this panel
    N = numel(dL)+1;

    %singularity treatment
    for i = 2:numel(dL)-1
        
       %distance form the singularity
       r1 = sqrt((x(1+6*(i-1):6*i,i)-x0(i)).^2+(y(1+6*(i-1):6*i,i)-y0(i)).^2);
       r2 = sqrt((x(1+6*(i-1):6*i,i)-x0(i+1)).^2+(y(1+6*(i-1):6*i,i)-y0(i+1)).^2);
       
       %single layer
       SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + 2*log(r1)-1;
       SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + 2*log(r2)-1;
       SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + 2*log(r1)-1;
       SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + 2*log(r2)-1;
       
    end
    
    %because I don't have to treat the singularity ON THE AXIS
    r = sqrt((x(1:6,2)-x0(2)).^2+(y(1:6,2)-y0(2)).^2);
    SXX(1:6,2) = SXX(1:6,2) + 2*log(r)-1;
    SYY(1:6,2) = SYY(1:6,2) + 2*log(r)-1;
    
    r = sqrt((x(1+6*(N-2):6*(N-1),N-1)-x0(N-1)).^2+(y(1+6*(N-2):6*(N-1),N-1)-y0(N-1)).^2);
    SXX(1+6*(N-2):6*(N-1),N-1) = SXX(1+6*(N-2):6*(N-1),N-1) + 2*log(r)-1;
    SYY(1+6*(N-2):6*(N-1),N-1) = SYY(1+6*(N-2):6*(N-1),N-1) + 2*log(r)-1;
    
    %analytical integration
    anA = diag(-dL(2:N-2).*log(dL(2:N-2))+2*dL(2:N-2));
    anB = diag(-dL(2:N-2).*log(dL(2:N-2))+dL(2:N-2));

end

