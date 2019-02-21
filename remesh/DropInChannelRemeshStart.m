%apply drop wall distribution for droplet in channel

function [a,b] = DropInChannelRemeshStart(a,b,PARAM)

    n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;
    walls = n+m+j;

    for i = 1:10
        
          %APPLY DISTRIBUTION ON DROPS AND WALL
          
          %compute distance between droplet 1 to wall
          a_this = a(walls+2:walls+2+q);
          b_this = b(walls+2:walls+2+q);
          a_other = a(n+1:n+m+1);
          b_other = b(n+1:n+m+1);

          Athis = repmat(a_this,numel(a_other),1);
          Bthis = repmat(b_this,numel(b_other),1);
          Aother = repmat(a_other,numel(a_this),1)';
          Bother = repmat(b_other,numel(b_this),1)';

          DIST = sqrt((Athis-Aother).^2 + (Bthis-Bother).^2);
          dist1 = min(DIST);
          dist3 = min(DIST');
          
          %apply distr to drop 1
          a_drop1 = a(walls+2:walls+2+q);
          b_drop1 = b(walls+2:walls+2+q);
          ratio = 1;

          [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric(a_drop1, b_drop1);
          K1 = curv_spline2(bx,by,cx,cy,dx,dy);
          N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2);...
                        -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2)];

          K2 = N(2,:)./b_drop1;
          K2([1 end]) = K1([1 end]);

          K = K1+K2;
          
          if PARAM.opt~=1

              distr1 = choose_distribution(a_drop1,b_drop1,K,PARAM.tune,q,dist1,PARAM.distr,ratio);
              [a(walls+2:walls+2+q), b(walls+2:walls+2+q)] = remesh_distribution(a_drop1,b_drop1,distr1);
          
          end
          
          %apply distr on wall
          a_wall = a(n+1:m+n+1);
          b_wall = b(n+1:m+n+1);
          
          distr_wall = choose_distribution(a_wall,b_wall,1,PARAM.tune_wall,PARAM.m,dist3,PARAM.distr,ratio);
          [a(n+1:n+m+1), b(n+1:n+m+1)] = remesh_distribution(a_wall,b_wall,distr_wall);
          
    end 

end