%nearly singular tretament of single layer potential with subtarction of
%curvature of the projected point as in Zichencko & Davis 2006

function z = nearly_sing_1layer(a,b,K,N,crit,PARAM,ax,ay,bx,by,cx,cy,dx,dy)

    K = repmat(K,numel(a),1);
    N1 = repmat(N(1,:),numel(a),1);
    N2 = repmat(N(2,:),numel(a),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% DROP 1 ON 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %AS A FIRST ATTEMPT I WILL CHOOSE ARBITRARLY THE POINTS WHERE TO DO
    %NEARLY SINGULAR TREATMENT
    ax2 = ax(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    ay2 = ay(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    bx2 = bx(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    by2 = by(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    cx2 = cx(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    cy2 = cy(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    dx2 = dx(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    dy2 = dy(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
       
    %compute splines coordinates
    t = 0:0.1:0.9;
    ttt = repmat(t,1,numel(ax2));
    axxx = reshape(repmat(ax2,numel(t),1),1,numel(ax2)*numel(t));
    bxxx = reshape(repmat(bx2,numel(t),1),1,numel(bx2)*numel(t));
    cxxx = reshape(repmat(cx2,numel(t),1),1,numel(cx2)*numel(t));
    dxxx = reshape(repmat(dx2,numel(t),1),1,numel(dx2)*numel(t));
    ayyy = reshape(repmat(ay2,numel(t),1),1,numel(ay2)*numel(t));
    byyy = reshape(repmat(by2,numel(t),1),1,numel(by2)*numel(t));
    cyyy = reshape(repmat(cy2,numel(t),1),1,numel(cy2)*numel(t));
    dyyy = reshape(repmat(dy2,numel(t),1),1,numel(dy2)*numel(t));
        
    %splines coordinates
    x = axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3;
    y = ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3;
    
    %compute the angular coefficient of the line perpendicular to the spline
    no = sqrt((bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2).^2+(byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2).^2);
    nx = (byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2)./no;
    ny = -(bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2)./no;
    m = ny./nx;
    
    %find out closest points (drop 1 on 2)  SIMPLER IF q=p
%     for i = round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1))
%         
%         %display(i)
%         
%         %distance from the node to the closest element of the other droplet
%         %[dist,ind] = min(sqrt((a(i)-a(q+2:end)).^2+(b(i)-b(q+2:end)).^2));
%         
%         %if dist<crit
%         
%         RES = b(i)-y-m.*(a(i)-x);
%         [MIN,ind] = min(abs(RES));
%         
% %         if MIN<10e-3
% %             warning('high residuals for the prjection point')
% %         end
%         
%         %projection point
%         %Xs = x(ind);
%         %Ys = y(ind);
%         K1s = ((bxxx(ind)+2*cxxx(ind)*ttt(ind)+3*dxxx(ind)*ttt(ind)^2)*(2*cyyy(ind)+6*dyyy(ind)*ttt(ind))-...
%             (byyy(ind)+2*cyyy(ind)*ttt(ind)+3*dyyy(ind)*ttt(ind)^2)*(2*cxxx(ind)+6*dxxx(ind)*ttt(ind)))/no(ind)^3;
%         K2s = ny(ind)/y(ind);
%         if isnan(K2s)==1
%             K2s = K1s;
%         end
%         Ks = K1s+K2s;
%         
%         K(i,PARAM.q+2:end) = K(i,PARAM.q+2:end) - Ks;
%         
% %         figure(10)
% %         plot(x,y)
% %         axis equal
% %         xlabel('x')
% %         ylabel('r')
% %         hold on
% %         plot(a(i),b(i),'or')
% %         plot([a(i) Xs],[b(i) Ys],'k')
% %         hold off
%         
%     end
    
    %copy for vectorial operation
    xxx = repmat(x',1,numel(a(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)))));
    yyy = repmat(y',1,numel(a(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)))));
    aaa = repmat(a(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1))),numel(x),1);
    bbb = repmat(b(round(numel(b(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1))),numel(x),1);
    mmm = repmat(m',1,numel(a(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)))));
    %find the line paerdendicular to the spline curve and passing for the
    %node
    RES = bbb-yyy-mmm.*(aaa-xxx);
    [~,ind] = min(abs(RES));
    ind(2:end) = ind(2:end) + cumsum(ones(1,numel(ind)-1));
    K1s = ((bxxx(ind)+2*cxxx(ind).*ttt(ind)+3*dxxx(ind).*ttt(ind).^2).*(2*cyyy(ind)+6*dyyy(ind).*ttt(ind))-...
            (byyy(ind)+2*cyyy(ind).*ttt(ind)+3*dyyy(ind).*ttt(ind).^2).*(2*cxxx(ind)+6*dxxx(ind).*ttt(ind)))./no(ind).^3;
    K2s = ny(ind)./y(ind);
    Ks = K1s+K2s;
    
    K(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)),PARAM.q+2:end) =...
        K(round(numel(a(1:PARAM.q+1))-PARAM.q*0.5):numel(a(1:PARAM.q+1)),PARAM.q+2:end) - ...
        repmat(Ks',1,PARAM.p+1);
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% DROP 2 ON 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %AS A FIRST ATTEMPT I WILL CHOOSE ARBITRARLY THE POINTS WHERE TO DO
    %NEARLY SINGULAR TREATMENT
    ax1 = ax(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    ay1 = ay(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    bx1 = bx(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    by1 = by(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    cx1 = cx(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    cy1 = cy(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    dx1 = dx(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    dy1 = dy(PARAM.q-round(PARAM.q*0.4):PARAM.q);
       
    %compute splines coordinates
    ttt = repmat(t,1,numel(ax1));
    axxx = reshape(repmat(ax1,numel(t),1),1,numel(ax1)*numel(t));
    bxxx = reshape(repmat(bx1,numel(t),1),1,numel(bx1)*numel(t));
    cxxx = reshape(repmat(cx1,numel(t),1),1,numel(cx1)*numel(t));
    dxxx = reshape(repmat(dx1,numel(t),1),1,numel(dx1)*numel(t));
    ayyy = reshape(repmat(ay1,numel(t),1),1,numel(ay1)*numel(t));
    byyy = reshape(repmat(by1,numel(t),1),1,numel(by1)*numel(t));
    cyyy = reshape(repmat(cy1,numel(t),1),1,numel(cy1)*numel(t));
    dyyy = reshape(repmat(dy1,numel(t),1),1,numel(dy1)*numel(t));
        
    %splines coordinates
    x = axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3;
    y = ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3;
    
    %compute the angular coefficient of the line perpendicular to the spline
    no = sqrt((bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2).^2+(byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2).^2);
    nx = (byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2)./no;
    ny = -(bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2)./no;
    m = ny./nx;
    
    %find out closest points (drop 1 on 2)  SIMPLER IF q=p
    for i = PARAM.q+2:PARAM.q+2+PARAM.p*0.5
        
        %display(i)
        
        %distance from the node to the closest element of the other droplet
        %[dist,ind] = min(sqrt((a(i)-a(q+2:end)).^2+(b(i)-b(q+2:end)).^2));
        
        %if dist<crit
        
        RES = b(i)-y-m.*(a(i)-x);
        [MIN,ind] = min(abs(RES));
        
%         if MIN<10e-3
%             warning('high residuals for the prjection point')
%         end
        
        %projection point
        %Xs = x(ind);
        %Ys = y(ind);
        K1s = ((bxxx(ind)+2*cxxx(ind)*ttt(ind)+3*dxxx(ind)*ttt(ind)^2)*(2*cyyy(ind)+6*dyyy(ind)*ttt(ind))-...
            (byyy(ind)+2*cyyy(ind)*ttt(ind)+3*dyyy(ind)*ttt(ind)^2)*(2*cxxx(ind)+6*dxxx(ind)*ttt(ind)))/no(ind)^3;
        K2s = ny(ind)/y(ind)*(ind>1)+K1s*(ind==1);
        if isnan(K2s)==1
            K2s = K1s;
        end
        Ks = K1s+K2s;
        
        K(i,1:PARAM.q+1) = K(i,1:PARAM.q+1) - Ks;
        
%         figure(10)
%         plot(x,y)
%         axis equal
%         xlabel('x')
%         ylabel('r')
%         hold on
%         plot(a(i),b(i),'or')
%         plot([a(i) Xs],[b(i) Ys],'k')
%         hold off
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %modified stresses
    df = [K(:,1:1+PARAM.q)*PARAM.gamma1 K(:,2+PARAM.q:end)*PARAM.gamma2];
    df_x = df.*N1;
    df_y = df.*N2;
    
    %matrix for BC
    z = zeros(2*(PARAM.q+PARAM.p+2));
    z(1:2:end-1,1:2:end-1) = df_x;
    z(2:2:end,1:2:end-1) = df_x;
    z(1:2:end-1,2:2:end) = df_y;
    z(2:2:end,2:2:end) = df_y;

end