%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [SXX,SYY,anA1x,anB1x,anA2x,anB2x,anA1y,anB1y,anA2y,anB2y] = stokesSTlinearSplineConvergence(SXX,SYY,eta,beta,hx,hy,h0x,h0y,h1x,h1y,x0,y0,GPt)

    % number of singularities on this panel
    N = numel(h0x)+1;

    %singularity treatment
    for i = 2:numel(h0x)-1
        
       %distance form the singularity
       r1 = sqrt((eta(1+6*(i-1):6*i)-x0(i)).^2+(beta(1+6*(i-1):6*i)-y0(i)).^2);
       r2 = sqrt((eta(1+6*(i-1):6*i)-x0(i+1)).^2+(beta(1+6*(i-1):6*i)-y0(i+1)).^2);
       
       %single layer
       SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + (2*log(r1).*hx(1+6*(i-1):6*i)...
           - 2*(log(r1./GPt).*hx(1+6*(i-1):6*i) + log(GPt).*(hx(1+6*(i-1):6*i)-h0x(i))))';
       SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + (2*log(r2).*hx(1+6*(i-1):6*i)...
           - 2*(log(r2./(1-GPt)).*hx(1+6*(i-1):6*i) + log(1-GPt).*(hx(1+6*(i-1):6*i)-h1x(i))))';
       
       SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + (2*log(r1).*hy(1+6*(i-1):6*i)...
           - 2*(log(r1./GPt).*hy(1+6*(i-1):6*i) + log(GPt).*(hy(1+6*(i-1):6*i)-h0y(i))))';
       SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + (2*log(r2).*hy(1+6*(i-1):6*i)...
           - 2*(log(r2./(1-GPt)).*hy(1+6*(i-1):6*i) + log(1-GPt).*(hy(1+6*(i-1):6*i)-h1y(i))))';
       
    end
    
    %because I don't have to treat the singularity ON THE AXIS
    r = sqrt((eta(1:6)-x0(2)).^2+(beta(1:6)-y0(2)).^2);
    SXX(1:6,2) = SXX(1:6,2) + (2*log(r).*hx(1:6) - 2*(log(r./(1-GPt)).*hx(1:6) + log(1-GPt).*(hx(1:6)-h1x(1))))';
    SYY(1:6,2) = SYY(1:6,2) + (2*log(r).*hy(1:6) - 2*(log(r./(1-GPt)).*hy(1:6) + log(1-GPt).*(hy(1:6)-h1y(1))))';
       
    r = sqrt((eta(1+6*(N-2):6*(N-1))-x0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-y0(N-1)).^2);
    SXX(1+6*(N-2):6*(N-1),N-1) = SXX(1+6*(N-2):6*(N-1),N-1) + (2*log(r).*hx(1+6*(N-2):6*(N-1)) - 2*(log(r./GPt).*hx(1+6*(N-2):6*(N-1)) + log(GPt).*(hx(1+6*(N-2):6*(N-1))-h0x(N-1))))';
    SYY(1+6*(N-2):6*(N-1),N-1) = SYY(1+6*(N-2):6*(N-1),N-1) + (2*log(r).*hy(1+6*(N-2):6*(N-1)) - 2*(log(r./GPt).*hy(1+6*(N-2):6*(N-1)) + log(GPt).*(hy(1+6*(N-2):6*(N-1))-h0y(N-1))))';
    
    %analytical integration
    anA1x = diag(3/2*h0x(2:N-2));
    anB1x = diag(1/2*h0x(2:N-2));
    anA2x = diag(3/2*h1x(2:N-2));
    anB2x = diag(1/2*h1x(2:N-2));
    anA1y = diag(3/2*h0y(2:N-2));
    anB1y = diag(1/2*h0y(2:N-2));
    anA2y = diag(3/2*h1y(2:N-2));
    anB2y = diag(1/2*h1y(2:N-2));

end

