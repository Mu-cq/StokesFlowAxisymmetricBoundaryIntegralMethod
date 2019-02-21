%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y,eta,beta]=drawline_adapt(x1,y1,x2,y2,N)

    %x length and y length of the line
    m=x2-x1;
    n=y2-y1;
    
    eta = linspace(x1,x2,N+1);
    beta = linspace(y1,y2,N+1);

    X = (1+tan(9*pi/16*(eta-x1)/m-9*pi/32)/tan(9*pi/32))*m/2+x1;
    Y = (1+tan(9*pi/16*(beta-y1)/n-9*pi/32)/tan(9*pi/32))*n/2+y1;
    
    if x1==x2
        X = x1*ones(numel(eta),1);
    end
    
    if y1==y2
        Y = y1*ones(numel(beta),1);
    end

%     figure
%     plot(X,Y,'o-')
%     axis equal
%     grid on
    
end