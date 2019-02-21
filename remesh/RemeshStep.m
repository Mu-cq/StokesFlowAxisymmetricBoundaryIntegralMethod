%add point s on wall when the interface comes very close, and remove when
%far

function [aNew,bNew,mNew] = RemeshStep(a,b,dist)

    m = numel(a)-1;

    %set when too distant
    FarCoeff = 2;
    %CloseCoeff = 0.5;

    %size of elements
    dl = sqrt(diff(a).^2+diff(b).^2);
    
    remesh = 1;
    %inefficient!!
    for i = 1:numel(dl)
        
        %remesh = i;
        
        %remesh when dist is too BIG compared to element size
        if dl(i) > dist(i)*FarCoeff
            
            m = m+1;
            %build new node
            a = [a(1:remesh) (a(remesh)+a(remesh+1))/2 a(remesh+1:end)];
            b = [b(1:remesh) (b(remesh)+b(remesh+1))/2 b(remesh+1:end)];
            
            remesh = remesh+1;
            
        end
        
        remesh = remesh+1;
        
    end
    
    aNew = a;   bNew = b; mNew = m;
    
end