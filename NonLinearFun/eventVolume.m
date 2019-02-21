%check if volume error is big and eventually repair

function [value,isterminal,direction] = eventVolume(t,xyDrop,V0,MaxErr)

    %drop coordinates
    x = xyDrop(1:2:end-1);  y = xyDrop(2:2:end);
    
    %drop volume
    V = axis_int_gauss_vect(x',y');
    
    %error on volume
    errV = abs(V-V0)/V0;
    
    value = 1;
    isterminal = 0;
    if errV>MaxErr
        
        %display('Error on volume')
        
        value = 0;
        isterminal = 1;
         
    end
        
    direction = 0;

end