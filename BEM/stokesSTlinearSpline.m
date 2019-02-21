%singular treatment for single and double layer potential for 2D Stokes: it
%gives the green function without the diverging part and analytically
%integrated part do add after the numerical integration

function [SXX,SYY,anA1,anB1,anA2,anB2] = stokesSTlinearSpline(SXX,SYY,eta,beta,h,h0,h1,x0,y0,GPt)

    % number of singularities on this panel
    N = numel(h0)+1;

    %singularity treatment
    for i = 2:numel(h0)-1
        
       %distance form the singularity
       r1 = sqrt((eta(1+6*(i-1):6*i)-x0(i)).^2+(beta(1+6*(i-1):6*i)-y0(i)).^2);
       r2 = sqrt((eta(1+6*(i-1):6*i)-x0(i+1)).^2+(beta(1+6*(i-1):6*i)-y0(i+1)).^2);
       
       %single layer
       SXX(1+6*(i-1):6*i,i) = SXX(1+6*(i-1):6*i,i) + (2*log(r1).*h(1+6*(i-1):6*i)...
           - 2*(log(r1./GPt).*h(1+6*(i-1):6*i) + log(GPt).*(h(1+6*(i-1):6*i)-h0(i))))';
       SXX(1+6*(i-1):6*i,i+1) = SXX(1+6*(i-1):6*i,i+1) + (2*log(r2).*h(1+6*(i-1):6*i)...
           - 2*(log(r2./(1-GPt)).*h(1+6*(i-1):6*i) + log(1-GPt).*(h(1+6*(i-1):6*i)-h1(i))))';
       
       SYY(1+6*(i-1):6*i,i) = SYY(1+6*(i-1):6*i,i) + (2*log(r1).*h(1+6*(i-1):6*i)...
           - 2*(log(r1./GPt).*h(1+6*(i-1):6*i) + log(GPt).*(h(1+6*(i-1):6*i)-h0(i))))';
       SYY(1+6*(i-1):6*i,i+1) = SYY(1+6*(i-1):6*i,i+1) + (2*log(r2).*h(1+6*(i-1):6*i)...
           - 2*(log(r2./(1-GPt)).*h(1+6*(i-1):6*i) + log(1-GPt).*(h(1+6*(i-1):6*i)-h1(i))))';
       
    end
    
    %because I don't have to treat the singularity ON THE AXIS
    r = sqrt((eta(1:6)-x0(2)).^2+(beta(1:6)-y0(2)).^2);
    SXX(1:6,2) = SXX(1:6,2) + (2*log(r).*h(1:6) - 2*(log(r./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1))))';
    SYY(1:6,2) = SYY(1:6,2) + (2*log(r).*h(1:6) - 2*(log(r./(1-GPt)).*h(1:6) + log(1-GPt).*(h(1:6)-h1(1))))';
       
    r = sqrt((eta(1+6*(N-2):6*(N-1))-x0(N-1)).^2+(beta(1+6*(N-2):6*(N-1))-y0(N-1)).^2);
    SXX(1+6*(N-2):6*(N-1),N-1) = SXX(1+6*(N-2):6*(N-1),N-1) + (2*log(r).*h(1+6*(N-2):6*(N-1)) - 2*(log(r./GPt).*h(1+6*(N-2):6*(N-1)) + log(GPt).*(h(1+6*(N-2):6*(N-1))-h0(N-1))))';
    SYY(1+6*(N-2):6*(N-1),N-1) = SYY(1+6*(N-2):6*(N-1),N-1) + (2*log(r).*h(1+6*(N-2):6*(N-1)) - 2*(log(r./GPt).*h(1+6*(N-2):6*(N-1)) + log(GPt).*(h(1+6*(N-2):6*(N-1))-h0(N-1))))';
    
    %analytical integration
    anA1 = diag(3/2*h0(2:N-2));
    anB1 = diag(1/2*h0(2:N-2));
    anA2 = diag(3/2*h1(2:N-2));
    anB2 = diag(1/2*h1(2:N-2));

end

