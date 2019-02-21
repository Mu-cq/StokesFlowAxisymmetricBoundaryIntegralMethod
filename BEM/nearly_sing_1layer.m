%nearly singular tretament of single layer potential with subtarction of
%curvature of the projected point as in Zichencko & Davis 2006

function z = nearly_sing_1layer(a,b,K,N,perc,PARAM,ax,ay,bx,by,cx,cy,dx,dy)

    %disp('ciao')

    %for vectorial operation
    K = repmat(K,numel(a),1);
    N1 = repmat(N(1,:),numel(a),1);
    %N2 = repmat(N(2,:),numel(a),1);
    N = N';
    N2 = repmat(N(:,2)',numel(a),1);
    
    %percentege of nodes on which I perfoem nearly singular treatment
    %perc = 0.2;
    
    %AS A FIRST ATTEMPT I WILL CHOOSE ARBITRARLY THE POINTS WHERE TO DO
    %NEARLY SINGULAR TREATMENT
    %choose the nodes on which performing the treatment
    range1 = round(numel(a(1:PARAM.q+1))-PARAM.q*perc):numel(a(1:PARAM.q+1));   %drop1
    range2 = PARAM.q+2:PARAM.q+2+round(PARAM.p*perc);   %drop2
    
    %spline coeff for drop1
    ax1 = ax(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    ay1 = ay(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    bx1 = bx(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    by1 = by(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    cx1 = cx(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    cy1 = cy(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    dx1 = dx(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    dy1 = dy(PARAM.q-round(PARAM.q*0.4):PARAM.q);
    
    %spline coeff for drop2
    ax2 = ax(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    ay2 = ay(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    bx2 = bx(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    by2 = by(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    cx2 = cx(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    cy2 = cy(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    dx2 = dx(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    dy2 = dy(PARAM.q+1:PARAM.q+round(PARAM.p*0.4));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% DROP 1 ON 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    x2 = axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3;
    y2 = ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3;
    
    %compute the angular coefficient of the line perpendicular to the spline
    no = sqrt((bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2).^2+(byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2).^2);
    nx = (byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2)./no;
    ny = -(bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2)./no;
    m = ny./nx;
    
    %here I fix the places where the nearly singular tretament is performed
    %copy for vectorial operation
    xxx = repmat(x2',1,numel(a(range1)));
    yyy = repmat(y2',1,numel(a(range1)));
    aaa = repmat(a(range1),numel(x2),1);
    bbb = repmat(b(range1),numel(x2),1);
    mmm = repmat(m',1,numel(a(range1)));
    %find the line paerdendicular to the spline curve and passing for the
    %node
    RES2 = bbb-yyy-mmm.*(aaa-xxx);
    [~,ind] = min(abs(RES2));
    ind(2:end) = ind(2:end) + cumsum(ones(1,numel(ind)-1));
    K1s = ((bxxx(ind)+2*cxxx(ind).*ttt(ind)+3*dxxx(ind).*ttt(ind).^2).*(2*cyyy(ind)+6*dyyy(ind).*ttt(ind))-...
            (byyy(ind)+2*cyyy(ind).*ttt(ind)+3*dyyy(ind).*ttt(ind).^2).*(2*cxxx(ind)+6*dxxx(ind).*ttt(ind)))./no(ind).^3;
    K2s = ny(ind)./y2(ind);
    Ks = K1s+K2s;
    
    K(range1,PARAM.q+2:end) = K(range1,PARAM.q+2:end) - repmat(Ks',1,PARAM.p+1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% DROP 2 ON 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
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
    x1 = axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3;
    y1 = ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3;
    
    %compute the angular coefficient of the line perpendicular to the spline
    no = sqrt((bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2).^2+(byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2).^2);
    nx = (byyy+2*cyyy.*ttt+3*dyyy.*ttt.^2)./no;
    ny = -(bxxx+2*cxxx.*ttt+3*dxxx.*ttt.^2)./no;
    m = ny./nx;
    
    %here I fix the places where the nearly singular tretament is performed
    %copy for vectorial operation
    xxx = repmat(x1',1,numel(a(range2)));
    yyy = repmat(y1',1,numel(a(range2)));
    aaa = repmat(a(range2),numel(x1),1);
    bbb = repmat(b(range2),numel(x1),1);
    mmm = repmat(m',1,numel(a(range2)));
    %find the line paerdendicular to the spline curve and passing for the
    %node
    RES1 = bbb-yyy-mmm.*(aaa-xxx);
    [~,ind] = min(abs(RES1));
    ind(2:end) = ind(2:end) + cumsum(ones(1,numel(ind)-1));
    K1s = ((bxxx(ind)+2*cxxx(ind).*ttt(ind)+3*dxxx(ind).*ttt(ind).^2).*(2*cyyy(ind)+6*dyyy(ind).*ttt(ind))-...
            (byyy(ind)+2*cyyy(ind).*ttt(ind)+3*dyyy(ind).*ttt(ind).^2).*(2*cxxx(ind)+6*dxxx(ind).*ttt(ind)))./no(ind).^3;
    K2s = ny(ind)./y1(ind);
    Ks = K1s+K2s;
    
    K(range2,1:PARAM.q+1) = K(range2,1:PARAM.q+1) - repmat(Ks',1,PARAM.q+1);
    
    %check perpendiculararity
%     for i = range1
%         
%         [MIN,ind] = min(abs(RES2));
%         Xs = x2(ind(i-range1(1)+1));
%         Ys = y2(ind(i-range1(1)+1));
%         
%         figure(10)
%         plot(x1,y1)
%         axis equal
%         xlabel('x')
%         ylabel('r')
%         hold on
%         plot(x2,y2)
%         plot(a(i),b(i),'or')
%         plot([a(i) Xs],[b(i) Ys],'k')
%         hold off
%         
%     end
%     
%     for i = range2
%         
%         [MIN,ind] = min(abs(RES1));
%         Xs = x1(ind(i-range2(1)+1));
%         Ys = y1(ind(i-range2(1)+1));
%         
%         figure(10)
%         plot(x1,y1)
%         axis equal
%         xlabel('x')
%         ylabel('r')
%         hold on
%         plot(x2,y2)
%         plot(a(i),b(i),'or')
%         plot([a(i) Xs],[b(i) Ys],'k')
%         hold off
%         
%     end
    
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