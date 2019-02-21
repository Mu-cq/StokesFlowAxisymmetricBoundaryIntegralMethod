function out = feval(f, xx, yy, zz)
%FEVAL   Evaluate a CHEBFUN3T object.
%    The workflow of this code is as follows:
%   1) Check the domain of f. If it is not the unit cube, then map the
%   points to the unit cube [-1, +1]^3. This should be done so that
%   Clenshaw's algorithm can be applied.
%
%   2) Check whether given points are vectors or tensors (Matrix case is 
%   rarely used and if we decided to add this at some point, we can copy 
%   from chebfun3/feval). If the input points are tensors generated by 
%   meshgrid or ndgrid, first determine the generating "vector" in each of 
%   them. Further operations should be applied only to these vectors and 
%   not to the full tensors (unless we have fully unstructured tensors of 
%   points).
% 
%   3) Call 1D Clenshaw's algorithm (usually only on the vectors because of 
%   Step 2).
%
%   4) Apply tensor-matrix contractions to get the output.
%
%   For a tensor of coefficients of size M*N*P and 3 vectors xeval, yeval, 
%   and zeval of size n*1, the complexity of this code is O (n*M*N*P).
%   Steps 1 and 2 can be reordered so that mapping is applied not to the
%   whole tensor, but to the generating vectors.
%
% See also CHEBFUN/FEVAL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

dom = f.domain;

%% Map all the points to [-1, +1]:
% If any of the 3 domains are not [-1, +1], then we need to map points in
% the corresponding domain [a, b] to [-1, +1] so that Clenshaw's algorithm 
% can be accurately applied. We are naiively following the work-flow 1D 
% Chebfun classes have.

% Map x coordinate of points to [-1, 1]:
if ( ~all(dom(1:2) == [-1, 1]) )
    % Create the linear map required in the 1st variable.
    linMap1 = bndfun.createMap(dom(1:2));
    % We just need the inverse part of the linMap1.
    xx = linMap1.Inv(xx);
end

% Map y coordinate of points to [-1, 1]:
if ( ~all(dom(3:4) == [-1, 1]) )
    % Create the linear map required in the 2nd variable.
    linMap2 = bndfun.createMap(dom(3:4));
    yy = linMap2.Inv(yy);
end

% Map z coordinate of points to [-1, 1]:
if ( ~all(dom(5:6) == [-1, 1]) )
    % Create the linear map required in the 3rd variable.
    linMap3 = bndfun.createMap(dom(5:6));
    zz = linMap3.Inv(zz);
end

%%
% %A question is how to apply Clenshaw's algorithm instead of the following 
% direct way ?
% This is done later in this code... Thanks to Anthony !
% Vmx = feval(chebpoly(0:N2-1),xeval);
% Vny = feval(chebpoly(0:M2-1),yeval);
% Vpz = feval(chebpoly(0:P2-1),zeval);
% simple_3D_eval = Vmx*(txm(simple_3D_coeffs,Vpz,3)*Vny.')
% % %simple_3D_eval = txm(txm(txm(simple_3D_coeffs,Vpz,3),Vny,2),Vmx,1)

n = size(xx);
[M,N,P] = size(f.coeffs);
if ( length(n) == 2 )
    % Inputs are vectors.
    %         Vmx = feval(chebpoly(0:M-1, dom(1:2)), xx);
    %         Vny = feval(chebpoly(0:N-1, dom(3:4)), yy);
    %         Vpz = feval(chebpoly(0:P-1, dom(5:6)), zz);
    % We should NOT use the above 3 lines because "chebpoly" explicitly
    % constructs several chebfun objects each time this is called. The 
    % following lines don't construct chebfuns at all and work just with 
    % vectors and matrices.
    % WARNING: This will NOT be very accurate for coordinate points way 
    % outside [-1,1], in which case the above lines provide more accuracy.
    Vmx = chebtech.clenshaw(xx, eye(M));
    Vny = chebtech.clenshaw(yy, eye(N));
    Vpz = chebtech.clenshaw(zz, eye(P));
    
    vals3D = zeros(n(1), 1);
    for i = 1:n(1) % number of rows in xx
        vals3D(i) = Vmx(i, :) * (chebfun3t.txm(f.coeffs, Vpz(i, :), 3) * ...
            Vny(i, :)');
    end
    
elseif ( length(n) == 3 )
    % Inputs are tensors.
    % If the evaluation points are derived from ndgrid, then there is a
    % fast way to evaluate a chebfun3. Check for this property. 
    if ( max(max(max(abs(bsxfun(@minus, xx, xx(:,1,1)))))) == 0  && ... 
            max(max(max(abs(bsxfun(@minus, yy, yy(1,:,1))))) == 0) && ... 
            max(max(max(abs(bsxfun(@minus, zz, zz(1,1,:))))) == 0) )
        xx = xx(:, 1, 1);
        yy = yy(1, :, 1).';
        zz = squeeze(zz(1, 1, :));
        Vmx = chebtech.clenshaw(xx, eye(M)); 
        Vny = chebtech.clenshaw(yy, eye(N));
        Vpz = chebtech.clenshaw(zz, eye(P));
        vals3D = chebfun3t.txm(chebfun3t.txm(chebfun3t.txm(f.coeffs, ...
            Vpz, 3), Vny, 2), Vmx, 1);
        
    elseif ( max(max(max(abs(bsxfun(@minus, xx, xx(1,:,1)))))) == 0  && ... 
            max(max(max(abs(bsxfun(@minus, yy, yy(:,1,1))))) == 0) && ... 
            max(max(max(abs(bsxfun(@minus, zz, zz(1,1,:))))) == 0) )
        
        xx = xx(1, :, 1).';
        yy = yy(:, 1, 1);
        zz = squeeze(zz(1, 1, :));
        
        Vmx = chebtech.clenshaw(xx, eye(M));
        Vny = chebtech.clenshaw(yy, eye(N));
        Vpz = chebtech.clenshaw(zz, eye(P));
        vals3D = chebfun3t.txm(chebfun3t.txm(chebfun3t.txm(f.coeffs, ...
            Vpz, 3), Vny, 2), Vmx, 1);
        vals3D = permute(vals3D, [2 1 3]);
        
    else
        % Inputs are not obtained by ndgrid or meshgrid.
        % Warning: This might be very slow!
        xx1 = chebfun3t.unfold(xx, 1)';
        yy1 = chebfun3t.unfold(yy, 1)';
        zz1 = chebfun3t.unfold(zz, 1)';
        fevalNew = zeros(n(2)*n(3), n(1));
        for i=1:n(1)
            fevalNew(:,i) = feval(f, xx1(:, i), yy1(:, i), zz1(:, i));
        end
        vals3D = reshape(fevalNew', n(1), n(2), n(3));
    end
end
out = vals3D;
end