%remesh in order to have more nodes toward the droplet tail

function [c, d, ds, int] = adapt(a,b,N)

    %total number of nodes
    %N = numel(a)-1;
    
    %display(N)
    
    ds = sqrt((a(1:end-1)-a(2:end)).^2+(b(1:end-1)-b(2:end)).^2);
    s = sum(ds);
    %theta = acos(ds./a(2:end));
    
    %display(sum(distr))
    
    int = zeros(1,numel(ds)+1);
    int(1) = 0;
    
    %compute the cumulative curve (integral of the distribution)
    for i=1:numel(ds)
        int(i+1) = sum(ds(1:i));
    end
    
%     if sum(distr)<N-10e-10*N || sum(distr)>N+10e-10*N
%         display('Warning: the number of nodes is not conserved in the distribution!');
%     end
     
     %c = zeros(1,numel(a));
     %d = zeros(1,numel(a));
     
     %decide for the spacing between the nodes (based on curvature value)
     
     space = (exp(0:s/N:s)-1)/(exp(s)-1)*s;
     
     c = interp1(int,a,space,'spline');
     %c = spline(int,(a(1:end-1)+a(2:end))/2,1:101);
     %d = ppval(pp,1:N);
     
     d = interp1(int,b,space,'spline');
     %d = spline(a,b,c);

end