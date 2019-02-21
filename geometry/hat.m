%given the physical domain, build the hat functions for linar elements
%already evaluated in the six Gauss point

function [phia, phib] = hat (l)

        %number of element
        N = numel(l);

        %Gauss points
        GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
        
        % points where I perform gauss integration
        GGP(1:6:6*N-5) = GP(1);
        GGP(2:6:6*N-4) = GP(2);
        GGP(3:6:6*N-3) = GP(3);
        GGP(4:6:6*N-2) = GP(4);
        GGP(5:6:6*N-1) = GP(5);
        GGP(6:6:6*N) = GP(6);

        %moltiplicate l
        ll(1:6:6*N-5) = l;
        ll(2:6:6*N-4) = l;
        ll(3:6:6*N-3) = l;
        ll(4:6:6*N-2) = l;
        ll(5:6:6*N-1) = l;
        ll(6:6:6*N) = l;
        
        %local frame of reference for every element
        s = (GGP+1)/2.*ll;
    
        %compute the values of the hat functions for the different elements
        phia = -1./ll.*s+1;
        phib = 1./ll.*s;
        
%         figure
%         plot(s,phia,'o-')
%         hold on
%         plot(s,phib,'or-')
%         hold off
    
return