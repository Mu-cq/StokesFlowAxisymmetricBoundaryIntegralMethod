%volume correction based on
%Pair collisions of fluid-filled elastic capsules in shear flow: Effects of membrane properties and polymer additives
%Pratik Pranay, Samartha G. Anekal, Juan P. Hernandez-Ortiz and Michael D. Graham 

function [a_out,b_out,A_after,P_after] = volume_correction2D(how_many,a_in,b_in,A_in,smooth)

        %compute area and volume
        P_before = perimeter2D(a_in,b_in);
        A_before = compute_area_2D(a_in',b_in');

        %%%%%%%%%%%%% VOLUME CORRECTION %%%%%%%%%%%%%%%%%%%
        
        for i = 1:how_many
            
            %interface normal and curvature
            [~, bx, ~, ~, ~, by, ~, ~] = my_spline_periodic (a_in, b_in);

            %compute the versor normal to the node N
            N = [by./sqrt(bx.*bx+by.*by);...
                  -bx./sqrt(bx.*bx+by.*by)];
              
            N = [N(1,:) N(1,1); N(2,:) N(2,1)];

            %display('Volume correction')
            a_in = a_in-(A_before-A_in)/P_before*N(1,:)*smooth;
            b_in = b_in-(A_before-A_in)/P_before*N(2,:)*smooth;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        A_before = compute_area_2D(a_in',b_in');
        P_before = perimeter2D(a_in,b_in);
        
        end
        
        A_after = A_before;
        P_after = P_before;
        a_out = a_in;
        b_out = b_in;

end

