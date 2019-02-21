% Solids_of_Revolution.m
%
% This script takes a function, rotates it about a specified axis, and calculates the volume of the solid
% The script is currently limited in that functions can only be entered in the form y=f(x). A new version will allow for entering of functions as x=f(y).
% Also you are limited to revolving around the x-axis (i.e. y=0) or lines parallel to the x-axis (i.e. y=k, where k can be poisitve or negative).
% A new version will allow for revolving around the y-axis (x=0) and lines parallel to it.
% However a current workaround is to switch all your x's and y's and perform the problem as such.  The shape will look the same but will be oriented
% differently.
%
% For example, to revolve the region bounded by y=x^2 and y=0 from x=0 to x=2 around the line x=5,
% revolve the region bounded by x=y^2 and y=2 around from x=0 to x=4 around the line y=5, but enter in x=y^2 as y=x^0.5.
% The volume, surface area, and arc lengths values will remain the same.
%
% Author: A. Bolu Ajiboye
% Created: Sunday, July 08, 2007
% Last Modified: Tuesday, July 10, 2007 by ABA



% Enter in function: Comment out either of the first five lines or the last five lines
% The first five lines prompt the user for the values from the command line
% the second five lines hard code the values, so no prompt is needed
my_function_out = input('Enter the outer function to revolve in the form y=f(x).\nFor the nth root, enter it in as "^(1/n) or as a decimal".  Example: y = 4*x^(3/2) - log(x^(1/2)) + e^(2*x) : ', 's');
my_function_in = input('\nEnter the inner function to revolve in the form y=f(x).\nFor the nth root, enter it in as "^(1/n) or as a decimal".  Example: y = 4*x^(3/2) - log(x^(1/2)) + e^(2*x) : ', 's');
min_x = input('\nWhat is the minimum x-value? ');
max_x = input('\nWhat is the maximum x-value? ');
my_axis = input('\nRevolve around (example: y = 5): ','s');
% my_function_out = 'y = x/4';
% my_function_in = 'y = x^0.5';
% min_x = 0;
% max_x = 16;
% my_axis = 'y = 5';

delta_x = (max_x-min_x)/75;
delta_y = delta_x;
delta_x_vol = delta_x/20000;      % If you get OUT OF MEMORY errors, then increase this value i.e. divide by a larger number

% Convert functions to an expression MATLAB can understand %
new_function_out = [];
for i = 1:length(my_function_out)
    if (my_function_out(i) == '*') || (my_function_out(i) == '/') || (my_function_out(i) == '^')
        new_function_out(end+1) = '.';
    end
    new_function_out(end+1) = my_function_out(i);
    if my_function_out(i) == 'e'
        new_function_out(end+1:end+5) = 'xp(1)';
    end
end
new_function_out = char(new_function_out);

new_function_in = [];
for i = 1:length(my_function_in)
    if (my_function_in(i) == '*') || (my_function_in(i) == '/') || (my_function_in(i) == '^')
        new_function_in(end+1) = '.';
    end
    new_function_in(end+1) = my_function_in(i);
    if my_function_in(i) == 'e'
        new_function_in(end+1:end+5) = 'xp(1)';
    end
end
new_function_in = char(new_function_in);
% End %

% Determine [x,y,z] points %
x = min_x:delta_x:max_x;
eval([new_function_out ';']);
if length(y) == 1
    y = y*ones(1,length(x));
end
y_out = abs(y - str2num(my_axis(find(my_axis == '=')+1:end)));
my_volume = [];
max_y = max(y_out);
max_z = -10000; min_z = 100000;
for i = 1:length(y_out)
    ytemp = y_out(i):-delta_y:-y_out(i);
    ztemp = sqrt(y_out(i)^2 - ytemp.^2);

    my_volume(i).x_out = x(i)*ones(1,2*length(ytemp));
    my_volume(i).y_out = [ytemp -ytemp];
    my_volume(i).z_out = [ztemp -ztemp]+str2num(my_axis(find(my_axis == '=')+1:end));
    max_z = max([max_z my_volume(i).z_out]);
    min_z = min([0 min_z my_volume(i).z_out]);
end

eval([new_function_in ';']);
if length(y) == 1
    y = y*ones(1,length(x));
end
y_in = abs(y - str2num(my_axis(find(my_axis == '=')+1:end)));
for i = 1:length(y_in)
    ytemp = y_in(i):-delta_y:-y_in(i);
    ztemp = sqrt(y_in(i)^2 - ytemp.^2);

    my_volume(i).x_in = x(i)*ones(1,2*length(ytemp));
    my_volume(i).y_in = [ytemp -ytemp];
    my_volume(i).z_in = [ztemp -ztemp]+str2num(my_axis(find(my_axis == '=')+1:end));
end
% End %

% Plot y = f(x) %
figure(1); set(gcf, 'windowstyle', 'docked', 'color', 'white');
eval([new_function_out ';']);
if length(y) == 1
    y = y*ones(1,length(x));
end
plot(x,y,'b-'); hold on;
eval([new_function_in ';']);
if length(y) == 1
    y = y*ones(1,length(x));
end
plot(x,y,'r-'); hold on;
max_y = max([max_y 0.000001]);
axis([min_x max_x -max_y max_y])
% legend([cellstr(my_function_out); cellstr(my_function_in)]);
% End %

% Plot volume by x-sections %
figure(2); set(gcf, 'windowstyle', 'docked', 'color', 'white');
for i = 1:length(y_out)
    plot3(my_volume(i).x_out,my_volume(i).y_out,my_volume(i).z_out,'b-'); hold on
end
for i = 1:length(y_in)
    plot3(my_volume(i).x_in,my_volume(i).y_in,my_volume(i).z_in,'r-'); hold on
end

xlabel('X','fontsize',14);
ylabel('Y','fontsize',14);
zlabel('Z','fontsize',14);

line([0 10], [0 0], [0 0], 'color', 'k', 'linewidth',2,'linestyle',':');
line([0 0], [0 10], [0 0], 'color', 'k', 'linewidth',2,'linestyle',':');
line([0 0], [0 0], [0 10], 'color', 'k', 'linewidth',2,'linestyle',':');
axis([min(min_x,0) max_x -max_y max_y min_z max_z]);
% End %

disp(' ');
disp(['Revolving the region bounded by ' my_function_out ' and ' my_function_in ' around ' my_axis ' from x = ' int2str(min_x) ' to x = ' int2str(max_x) '...']);
