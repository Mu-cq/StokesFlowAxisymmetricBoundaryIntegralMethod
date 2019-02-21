%volume correction by resizing

function [x,y] = volumeCorrectionResize(x,y,V0)

    %compute current volume
    V = axis_int_gauss_vect(x,y);
    
    %compute how much to resize
    alpha = nthroot(V0/V,3);
    
    %resize
    x = x*alpha;
    y = y*alpha;

end