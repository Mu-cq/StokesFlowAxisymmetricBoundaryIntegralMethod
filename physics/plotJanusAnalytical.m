%plot species concentration for a Janus Particle

function [c,DpsiDr] = plotJanusAnalytical(r,theta,thetaC)

%parameters
muC = cos(thetaC);
nTruncation = 500;                %number of polynomials

%mobility
M = -1;
D1 = finiteDifference1D(numel(theta),[2 0],1);

%compute concentration
k0 = (1-muC)/2;
c = -k0./r;     % firt term of the series
DpsiDr = 0;
for i = 1:nTruncation
    
    %compute Legendre Polynomials
    P = legendre(i+1,muC);
    Pnplus1 = P(1,:)';
    P = legendre(i-1,muC);
    Pnminus1 = P(1,:)';
    P = legendre(i,cos(theta));
    Pn = P(1,:)';
    PnDer = -D1*Pn./sin(theta')/pi;
    PnDer([1 end]) = i*(i+1)/2;
    if sum(i == 2:2:nTruncation)
        PnDer(end) = -i*(i+1)/2;
    end
    
    %add term to the serie
    kn = 0.5*(Pnminus1-Pnplus1);
    cp = -kn/(i+1)./r.^(i+1);
    c = c + cp.*Pn;
    
    %series for slip velocity
    n = i;
    if i==1
        DpsiDrN = -(2*r^3+1)/(3*r^2);
    else
        DpsiDrN = 0.5*(-n/r^(n+1) + (n-2)/(r^(n-1)));
    end
    alphaN = -n*(n+1)*M/(2*n+1)*cp(1);
    %alphaN = i*kn*M/(2*i+1);
    DpsiDr = DpsiDr + alphaN*DpsiDrN*sin(theta').*PnDer*(2*n+1)/(n*(n+1));
    
end