%compute F, the function to minimize in order not to have mesh distorsion

function F = computeFclosed(x,h2,V)

    vij = [V(2:end,:); V(1,:)]-V;
    
    H = 4*(1./h2-h2./(sum(x.*x,2).*sum(x.*x,2))).^2;
    
    f = H.*(sum(x.*vij,2)).^2;
    
    F = sum(f);

end