%convert Gerris data to my geometry data

function [Xdrop,Ydrop] = ConvertGerris(x,y)
    
    %delete duplicate
    [x,ind] = unique(x);    y = y(ind);
    %[y,ind] = unique(y);    x = x(ind);
    
    %identify negative numbers of y
    ind = find(y>=0);
    
    [Xpos,ind2] = sort(x(ind),'descend');
    
    yyy = y(ind);
    Ypos = yyy(ind2);
    
    %by hand
    Xpos = [Xpos(1); Xpos];
    Ypos = [0; Ypos];
    
%     figure
%     plot(Xpos,Ypos,'o-')
%     axis equal
%     grid on
    
    %mirror negative values
    Xdrop = [Xpos; flip(Xpos(1:end-1))];
    Ydrop = [Ypos; -flip(Ypos(1:end-1))];
    
    %resize for matlab
    xcm = center_mass_2D(Xdrop,Ydrop);
    Xdrop = Xdrop-xcm;
    Xdrop = 2*Xdrop;    Ydrop = 2*Ydrop;
    
%     figure
%     plot(Xdrop,Ydrop,'o-')
%     axis equal
%     grid on

end