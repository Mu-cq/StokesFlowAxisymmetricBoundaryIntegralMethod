%remesh using a chosen distribution, this also can change the number of
%elements: new number of elements is equal to the numbear od data points in
%the distribution

function distr = choose_distribution(a,b,K,tune,q,dist,opt,ratio)

    %average because I want a distribution on the spacing
    K = (K(1:end-1)+K(2:end))/2;
    dist = (dist(1:end-1)+dist(2:end))/2;
    Xsing = (a(1:end-1)+a(2:end))/2;
    Ysing = (b(1:end-1)+b(2:end))/2;

    if opt==1   %remesh based on curvature
        
        %choose distribution as the reciprocal of the curvature, the tune
        %parameter smooth the distribution when is larger (adding a contsnat value the difference between maximum and minimum is smaller)
        temp = 1./(abs(K)+max(abs(K))*tune)';
        
    elseif opt==2   %remesh based on distance from wall
        
        %this is in order to have the same weigth for the 2 distribution or
        %detrmined by ratio
        temp = (dist+max(dist)*tune)';
        
    elseif opt==3   %remesh based on curvature and distance from wall
        
        %this is in order to have the same weigth for the 2 distribution or
        %detrmined by ratio
        areaK = sum(1./K); areaDIST = sum(dist);
        dist = dist*(areaK/areaDIST)*ratio;
        temp = (dist+max(dist)*tune)' + 1./(abs(K)+max(abs(K))*tune)';
        
    elseif opt==4   %more mesh close to axis
        
        temp = Ysing'+max(Ysing)*tune;
        
    elseif opt==5   %more mesh on the left part of the droplet
        
        temp = Xsing'-min(a) + max(Xsing'-min(Xsing)*tune);
        
    elseif opt==6   %more mesh on the rigth part of the droplet
        
        temp = Xsing'-min(Xsing) + max(Xsing'-min(Xsing)*tune);
        temp = 1./temp;
        
    elseif opt==7   %more mesh on the rigth part of the droplet MORE THAN 6
        
        %distance from right part
        distRight = max(Xsing)-Xsing;
        
        %distribution
        temp = 1./(exp(distRight' + tune)-1);
        temp = 1./temp;
        
%     elseif opt==8   %dedicated to edge state (constant on tail and than smoooth)
%         
%         %distance form right part
%         distRight = max(Xsing)-Xsing;
%         
%         %distribution

    elseif opt==8     % uniform mesh on droplet
        
        temp = ones(numel(K),1);
        
    end
    
    %fit smoothing spline and evaluete data where needed
    f = fit((1:numel(temp))',temp,'smoothingspline');
    distr = feval(f,(1:q)'/q*numel(temp));

end