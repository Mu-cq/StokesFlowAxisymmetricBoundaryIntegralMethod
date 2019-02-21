%remesh using a chosen distribution, this also can change the number of
%elements: new number of elements is equal to the numbear od data points in
%the distribution

function distr = choose_distribution_nodes(a,b,K,tune,nodes,dist,opt,ratio)

    %xcm = center_mass(a,b);

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
        
        temp = b+max(b)*tune;
        
    elseif opt==5   %more mesh on the left part of the droplet
        
        temp = a-min(a) + max(a-min(a))*tune;
        temp = 1./temp;
        
    elseif opt==6   %more mesh on the rigth part of the droplet
        
        temp = a-min(a) + max(a-min(a))*tune;
        
    end
    
    figure(10)
    plot(temp)
    grid on
    
    %fit smoothing spline and evaluete data where needed
    f = fit((1:numel(temp))',temp,'smoothingspline');
    distr = feval(f,(1:nodes)'/nodes*numel(temp));

end