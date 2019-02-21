%volume correction based on
%Pair collisions of fluid-filled elastic capsules in shear flow: Effects of membrane properties and polymer additives
%Pratik Pranay, Samartha G. Anekal, Juan P. Hernandez-Ortiz and Michael D. Graham 

function [a_out,b_out,V_after,Area_after] = volume_correction(how_many,a_in,b_in,V_in,smooth)

        %compute area and volume
        Area_before = surf_gauss_vect(a_in,b_in);
        V_before = axis_int_gauss_vect(a_in,b_in);

        %%%%%%%%%%%%% VOLUME CORRECTION %%%%%%%%%%%%%%%%%%%
        
        for i = 1:how_many
            
            [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric (a_in, b_in);

            %compute the versor normal to the node N
            N = [by./sqrt(bx.*bx+by.*by) by(end)+2*cy(end)+3*dy(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
                  -bx./sqrt(bx.*bx+by.*by) -bx(end)-2*cx(end)-3*dx(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];

              %those because the points are on the axis
            N(:,1) = [1; 0];
            N(:,end) = [-1; 0];

            %display('Volume correction')
            a_in = a_in-(V_before-V_in)/Area_before*N(1,:)*smooth;
            b_in = b_in-(V_before-V_in)/Area_before*N(2,:)*smooth;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        V_before = axis_int_gauss_vect(a_in,b_in);
        Area_before = surf_gauss_vect(a_in,b_in);
        
        end
        
        V_after = V_before;
        Area_after = Area_before;
        a_out = a_in;
        b_out = b_in;

end

