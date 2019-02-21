%compute series of legendre form given data points

function x = LegendreFourierBuild2DgridInterp(u,v,fff,aaan,bbbn,nT,nS,PPP)

    x = zeros(numel(u),1);

    %compute x at (u,v)
    for i = 1:numel(u)
        
        %loop over Legendre
        for k = 1:nT
            
            P = PPP(k,i);
            
            %loop over Fourier
            for l = 1:nS
                
                %Legendre coefficient
                f = fff(k,l);
                an = aaan(k,l);
                bn = bbbn(k,l);
                
                if k==1
                    an = an/2;
                end
                
                x(i) = x(i) + f*P*(an*cos(l*v(i))+bn*sin(l*v(i)));
                
            end
            
        end
        
    end

end