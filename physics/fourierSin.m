%compute cosine modes coeffiecients

function f = fourierSin(theta,x,modes)

    f = zeros(modes,1);
    for i = 1:modes
        
        coeff = i-1;

        f(i) = 1/pi*trapz(theta,x.*sin(coeff*theta));
        
    end

end