function [F,E]= ell_int_vect (K)

%-----------------------------------------
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%--------------------------------------------
%  COMPUTATION OF COMPLETE ELLIPTIC INTEGRALS
%  OF THE FIRST AND SECOND KIND
%--------------------------------------------

tol=0.0000000000001;

%-----------
% iterations
%-----------

K2 = K.*K;

F = 0.5*pi;
E = 1.0;
G = 1.0;
B = K;
D = 10*tol;

%pause
%disp('one')

while max(max(abs(D))) > tol,
    %disp(max(max(D)))
  C = sqrt(1.0-B.^2);
  B = (1.0-C)./(1.0+C);
  D = F.*B;
  F = F+D;
  G = 0.50*G.*B;
  E = E+G;
%pause
%disp('two')
end

E = F.*(1.0-0.50*K2.*E);

%-----
% done
%-----

return
