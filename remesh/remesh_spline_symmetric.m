%add one nodes if element becomes too big copared to the curvature, viceversa
%one node less if element becomes too small

%BUG WHEN THE DROPLET IS FLOATING IN THE POSITIVE DIRECTION!!!

function [c, d, element] = remesh_spline_symmetric(a,b,min)

    %total number of elements
    N = numel(a)-1;
    element = N;
    
    %[normal,~]=normal_versor_clean(a,b,N);
    %radius = CurvDrop2([a(2) a a(end-1)], [-b(2) b -b(end-1)], [1 normal(1,2:end-1) -1; 0 normal(2,2:end-1) 0]);
  
    curvature = curv_spline(a,b);

    %first component of the curvature
    %K1 = (1./radius(1:end-1)+1./radius(2:end))/2;
    K1 = (curvature(1:end-1)+curvature(2:end))/2;
    

    %sum the two component of the curvature
    curv = K1;
    
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);
      
    
    ds = sqrt((a(1:end-1)-a(2:end)).^2+(b(1:end-1)-b(2:end)).^2);
    %s = sum(ds);
    c = zeros(1,10*N);
    d = zeros(1,10*N);
    count = 2;
    
    %control if add or remove elements
    control = 0;
    %active = 0;
    
    for i=1:N
        
        if abs(ds(i)*curv(i)) > 10*min %abs(ds(i)*curv(i)) > 15*min
            
            %Add element angle
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
            %active = 1;
            
        elseif ds(i) > 3*min %high
            
            %Add element absolute value
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
            %active = 1;

        elseif (ds(i) < min/1.3) %(ds(i) < min/1.1) 
            
            %Remove element
            element = element-1;
            %interp with spline

            x = ax(i)+bx(i)*0.5+cx(i)*0.5^2+dx(i)*0.5^3;
            y = ay(i)+by(i)*0.5+cy(i)*0.5^2+dy(i)*0.5^3;
            
            c(count-1) = x;
            d(count-1) = y;
            
            if i == N
                c(count-1) = ax(i)+bx(i)+cx(i)+dx(i);
                d(count-1) = ay(i)+by(i)+cy(i)+dy(i);
            end
            
            control = 1;
        else
            c(count) = a(i+1);
            d(count) = b(i+1);
            count = count+1;
            %active = 0;
        end

        
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