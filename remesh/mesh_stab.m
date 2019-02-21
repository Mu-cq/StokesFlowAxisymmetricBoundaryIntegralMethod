%passive mesh stabilization algorithm from 'Emulsion flow through a packed
%bed with,multiple drop breakup' Zinchenko

%for a 1D maesh every nodes is connected only with the one before and the
%one after, in my case i is equal to Nel-2 and j=2

function [V,check,ite_plus,saveF] = mesh_stab(x,y,v,k,k2,n,PARAM,distr)

    k2([1 end]) = k([1 end])/2;
    k([1 end]) = k([1 end])/2;
    
    %adaptivity parameters
    g = PARAM.g;
    ratio = PARAM.ratio;
    tune = PARAM.tune;
    dist_wall = PARAM.dist_wall;
    
    %choose option
    option = PARAM.opt;
    
    %number of edges
    N = numel(x)-1;
    
    %number of treated nodes
    tr = numel(x)-2;

    saveF = zeros(1000,1);
    
    %I apply the algorithm at the nodes not on the axis
    %the distance form the preceding node
    
    %%%%%%%%%%%%%%%% DEFINITION OF GEOMETRICAL VARIABLES %%%%%%%%%%%%%%%%
    
    xij = [x(2:end)-x(1:end-1) y(2:end)-y(1:end-1)];
    %CHECK THE COEFFICIENT HERE because I use only one component of the
    %curvature!!! 0.004 ok when the starting radius is 1!!!
    
    if option==1
        G = k.*k + k2.*k2 + 0.004; %THIS1
    elseif option==2
        G = k.*k + [1/y(2)^2 1./y(2:end-1)'.^2 1/y(end-1)^2]; %THIS2
    elseif option==3 %cluster on the left side
        x1 = min(x);
        %display('ciao1')
        cluster = 1./(x-x1+0.004).^10;
        G = cluster'; %THIS3
    elseif option==4 %cluster on the right side
        x1 = max(x);
        %display('ciao2')
        cluster = 1./(x1-x+0.004).^10;
        G = cluster'; %THIS3
    elseif option == 5
        %display('ciao')
        G = 1./choose_distribution_nodes(x,y,k+k2,tune,numel(k),dist_wall,distr,ratio)';
    %elseif option==6
     %   G = 1./PARAM.WG';
    end
    
    dl = sqrt(xij(:,1).^2+xij(:,2).^2);
    %this change respect to the paper because I'm in 2D
    K = 1/2/N*sum(G'.^g.*[dl(1).^2; dl(1:end-1).^2+dl(2:end).^2; dl(end).^2]);
    h2 = K*G'.^(-g);
    hij2 = 0.5*(h2(2:end)+h2(1:end-1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %check
    check = max(dot(xij,xij,2)./hij2)/min(dot(xij,xij,2)./hij2);
    
    %%%%%%%%%%%%%%%%%%% STEEPEST DESCENT METHOD %%%%%%%%%%%%%%%%%%%%%%%%
    
    vij = v(2:end,:)-v(1:end-1,:);
    F_old = computeF(xij,hij2,v);
     
    f = zeros(tr+2,2);
    
    H = 4*(1./hij2-hij2./(sum(xij.*xij,2).*sum(xij.*xij,2))).^2;
    
    f(2:end-1,:) = [2*H(1:end-1).*sum(xij(1:end-1,1).*vij(1:end-1,1),2).*xij(1:end-1,1) - 2*H(2:end).*sum(xij(2:end,1).*vij(2:end,1),2).*xij(2:end,1) ...
        2*H(1:end-1).*sum(xij(1:end-1,2).*vij(1:end-1,2),2).*xij(1:end-1,2) - 2*H(2:end).*sum(xij(2:end,2).*vij(2:end,2),2).*xij(2:end,2)];
    
    %take only the tangential velocities
    ft = zeros(tr+2,2);
    ft(:,1) = f(:,1)-sum(f.*n',2).*n(1,:)';
    ft(:,2) = f(:,2)-sum(f.*n',2).*n(2,:)';
    f = ft;
    
    fij = f(2:end,:)-f(1:end-1,:);
    
    %parabola coeffiecient
    a = sum(H .* sum(xij.*fij,2).^2);
    b = 2*sum(H .* sum(xij.*vij,2) .* sum(xij.*fij,2));
    
    xi = -b/2/a;
    
    V = v + xi*f;
    
    F = computeF(xij,hij2,V);
    
    if F>F_old
       warning('Object function F is not decreasing') 
    end
    
    saveF(1) = F_old;
    saveF(2) = F;
    
    v1 = v;
    v = V;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    w = 0;
    
    %%%%%%%%%%%%%%%%%%%%%% CONJUGATE GRADIENT METHOD %%%%%%%%%%%%%%%%%%%
    count = 0;
    
    while abs(F-F_old) > abs(F)*1e-6
        
        w = w+1;

        F_old = F;

        vij = v(2:end,:)-v(1:end-1,:);
        
        dv = v-v1;
        
        dvij = dv(2:end,:)-dv(1:end-1,:);

        f = zeros(tr+2,2);
        
        f(2:end-1,:) = [2*H(1:end-1).*sum(xij(1:end-1,1).*vij(1:end-1,1),2).*xij(1:end-1,1) - 2*H(2:end).*sum(xij(2:end,1).*vij(2:end,1),2).*xij(2:end,1) ...
            2*H(1:end-1).*sum(xij(1:end-1,2).*vij(1:end-1,2),2).*xij(1:end-1,2) - 2*H(2:end).*sum(xij(2:end,2).*vij(2:end,2),2).*xij(2:end,2)];

        %take only the tangential velocities
        ft = zeros(tr+2,2);
        ft(:,1) = f(:,1)-sum(f.*n',2).*n(1,:)';
        ft(:,2) = f(:,2)-sum(f.*n',2).*n(2,:)';
        f = ft;

        fij = f(2:end,:)-f(1:end-1,:);
        
        %paraboloid coeffiecient (as in Lailai notation)
        
        a = sum(H .* sum(xij.*fij,2).^2);
        c = sum(H .* sum(xij.*dvij,2).^2);
        b = 2*sum(H .* sum(xij.*fij,2) .* sum(xij.*dvij,2));
        d = 2*sum(H .* sum(xij.*vij,2) .* sum(xij.*fij,2));
        e = 2*sum(H .* sum(xij.*vij,2) .* sum(xij.*dvij,2));

        xi = (-2*c*d+b*e)/(4*a*c-b*b);
        eta = (-2*a*e+b*d)/(4*a*c-b*b);
        
        V = v + xi*f + eta*dv;
        
        F = computeF(xij,hij2,V);

        if F>F_old
            warning('Object function F is not decreasing') 
        end
        
%         if count>100
%             warning('Passive mesh stabilization didn't reach convergence')
%             break
%         end
        
        saveF(w+2) = F;
        
        %plot(F,'o')

        v1 = v;
        v = V;
        
        %display(count)
        
        count = count+1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ite_plus = w+2;
    saveF = saveF(1:w+2,1);
    
    %disp(count)

end

