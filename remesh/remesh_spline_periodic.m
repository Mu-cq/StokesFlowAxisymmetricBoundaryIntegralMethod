%add one nodes if element becomes too big copared to the curvature, viceversa
%one node less if element becomes too small

%BUG WHEN THE DROPLET IS FLOATING IN THE POSITIVE DIRECTION!!!

function [c, d, element] = remesh_spline_periodic(a,b,min)

    %total number of elements
    N = numel(a)-1;
    element = N;
  
    %interface normal and curvature
    [ax, bx, cx, dx, ay, by, cy, dy] = my_spline_periodic (a, b);
  
    %compute normals to the interface and curvature
    K = curv_spline_closed(bx,by,cx,cy);
    curv = ([K(1:end) K(1)]+[K(2:end) K(1) K(2)])/2;
      
    
    ds = sqrt((a(1:end-1)-a(2:end)).^2+(b(1:end-1)-b(2:end)).^2);
    %s = sum(ds);
    c = zeros(1,10*N);
    d = zeros(1,10*N);
    count = 2;
    
    %control if add or remove elements
    control = 0;
    active = 0;
    
    for i=1:N
        
        if abs(ds(i)*curv(i))>1
            
            %display('Add element angle')
            element = element+1;
            c(count+1) = a(i+1);
            d(count+1) = b(i+1);
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            x = ax(i)+bx(i)*0.5+cx(i)*0.5^2+dx(i)*0.5^3;
            y = ay(i)+by(i)*0.5+cy(i)*0.5^2+dy(i)*0.5^3;
            
            c(count) = x;
            d(count) = y;
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            count = count+2;
            control = 1;
            active = 1;
            
        elseif ds(i) > 5*min %ds(i) > 0.3 %high
            
            %display('Add element absolute value')
            element = element+1;
            c(count+1) = a(i+1);
            d(count+1) = b(i+1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            x = ax(i)+bx(i)*0.5+cx(i)*0.5^2+dx(i)*0.5^3;
            y = ay(i)+by(i)*0.5+cy(i)*0.5^2+dy(i)*0.5^3;
            
            c(count) = x;
            d(count) = y;
            
            %%%%%%%%%%%%%%%%%%%%%%%%
           
            count = count+2;
            control = 1;  
            active = 1;

        elseif (ds(i) < min/8)
            
            %display('Remove element')
            element = element-1;
            %interp with spline

            x = ax(i)+bx(i)*0.5+cx(i)*0.5^2+dx(i)*0.5^3;
            y = ay(i)+by(i)*0.5+cy(i)*0.5^2+dy(i)*0.5^3;
            
            c(count-1) = x;
            d(count-1) = y;
            
            control = 1;
        else
            c(count) = a(i+1);
            d(count) = b(i+1);
            count = count+1;
        end

        active = 0;
    end

    %set the first point
    c = c(1:element+1);
    d = d(1:element+1);
    c(1) = a(1);
    d(1) = b(1);

    if control
        display(strcat('Remesh: number of element is equal to ',num2str(element)))
    end

end