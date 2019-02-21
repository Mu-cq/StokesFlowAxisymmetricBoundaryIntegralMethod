%nearly singular tretament of single layer potential with subtraction of
%curvature of the projected point as in Zichencko & Davis 2006

function z = nearly_sing_1layer_WALL_DROP(X_wall,Y_wall,K,N,ax,ay,bx,by,cx,cy,dx,dy,gamma)

    K = repmat(K,numel(X_wall),1);
    N1 = repmat(N(1,:),numel(X_wall),1);
    N2 = repmat(N(2,:),numel(X_wall),1);
    
    %Kmod = zeros(numel(K));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% WALL ON DROPLET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    %compute splines coordinates
    t = 0:0.1:0.9;
    ttt = repmat(t,1,numel(ax));
    axxx = reshape(repmat(ax,numel(t),1),1,numel(ax)*numel(t));
    bxxx = reshape(repmat(bx,numel(t),1),1,numel(bx)*numel(t));
    cxxx = reshape(repmat(cx,numel(t),1),1,numel(cx)*numel(t));
    dxxx = reshape(repmat(dx,numel(t),1),1,numel(dx)*numel(t));
    ayyy = reshape(repmat(ay,numel(t),1),1,numel(ay)*numel(t));
    byyy = reshape(repmat(by,numel(t),1),1,numel(by)*numel(t));
    cyyy = reshape(repmat(cy,numel(t),1),1,numel(cy)*numel(t));
    dyyy = reshape(repmat(dy,numel(t),1),1,numel(dy)*numel(t));
        
    %splines coordinates
    x = axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3;
    y = ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3;
    
    %compute the angular coefficient of the line perpendicular to the spline
    no = sqrt((bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2).^2+(byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2).^2);
    nx = (byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2)./no;
    ny = -(bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2)./no;
    m = ny./nx;
    
    %find out closest points (betwen drop and wall)
    for i = 1:numel(X_wall)
        
        %display(i)
        
        RES = Y_wall(i)-y-m.*(X_wall(i)-x);
        [~,ind] = min(abs(RES));
        
        %projection point
        K1s = ((bxxx(ind)+2*cxxx(ind)*ttt(ind)+3*dxxx(ind)*ttt(ind)^2)*(2*cyyy(ind)+6*dyyy(ind)*ttt(ind))-...
            (byyy(ind)+2*cyyy(ind)*ttt(ind)+3*dyyy(ind)*ttt(ind)^2)*(2*cxxx(ind)+6*dxxx(ind)*ttt(ind)))/no(ind)^3;
        K2s = ny(ind)/y(ind);
        if isnan(K2s)==1
            K2s = K1s;
        end
        Ks = K1s+K2s;
        
        K(i,:) = K(i,:) - Ks;
        
%         figure(10)
%         plot(x,y)
%         axis equal
%         xlabel('x')
%         ylabel('r')
%         hold on
%         plot(X_wall(i),Y_wall(i),'or')
%         %plot([a(i) Xs],[b(i) Ys],'k')
%         hold off
        
    end
    
    %copy for vectorial operation
%     xxx = repmat(x',1,numel(a(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)))));
%     yyy = repmat(y',1,numel(a(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)))));
%     aaa = repmat(a(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1))),numel(x),1);
%     bbb = repmat(b(round(numel(b(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1))),numel(x),1);
%     mmm = repmat(m',1,numel(a(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)))));
%     %find the line paerdendicular to the spline curve and passing for the
%     %node
%     RES = bbb-yyy-mmm.*(aaa-xxx);
%     [~,ind] = min(abs(RES));
%     ind(2:end) = ind(2:end) + cumsum(ones(1,numel(ind)-1));
%     K1s = ((bxxx(ind)+2*cxxx(ind).*ttt(ind)+3*dxxx(ind).*ttt(ind).^2).*(2*cyyy(ind)+6*dyyy(ind).*ttt(ind))-...
%             (byyy(ind)+2*cyyy(ind).*ttt(ind)+3*dyyy(ind).*ttt(ind).^2).*(2*cxxx(ind)+6*dxxx(ind).*ttt(ind)))./no(ind).^3;
%     K2s = ny(ind)./y(ind);
%     Ks = K1s+K2s;
%     
%     K(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)),PARAM.q+2:end) =...
%         K(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)),PARAM.q+2:end) - ...
%         repmat(Ks',1,PARAM.p+1);
    
%     for i = round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1))
%         
%         [MIN,ind] = min(abs(RES));
%         Xs = x(ind(i-round(numel(a(1:PARAM.q+1))-PARAM.q*0.5)+1));
%         Ys = y(ind(i-round(numel(a(1:PARAM.q+1))-PARAM.q*0.5)+1));
%         
%         figure(10)
%         plot(x,y)
%         axis equal
%         xlabel('x')
%         ylabel('r')
%         hold on
%         plot(a(i),b(i),'or')
%         plot([a(i) Xs],[b(i) Ys],'k')
%         hold off
%         
%     end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %modified stresses
    df = K*gamma;
    df_x = df.*N1;
    df_y = df.*N2;
    
    %matrix for BC
    z = zeros(2*numel(X_wall),2*numel(ax)+2);
    z(1:2:end-1,1:2:end-1) = df_x;
    z(2:2:end,1:2:end-1) = df_x;
    z(1:2:end-1,2:2:end) = df_y;
    z(2:2:end,2:2:end) = df_y;

end