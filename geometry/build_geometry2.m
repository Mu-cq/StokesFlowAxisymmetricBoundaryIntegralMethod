%build geometry following the panel instructions

function PARAM = build_geometry2(PARAM)

    PARAM.varlistA = matlab.lang.makeUniqueStrings({'a','a','a','a','a','a','a','a','a'});
    PARAM.varlistB = matlab.lang.makeUniqueStrings({'b','b','b','b','b','b','b','b','b'});

    for i = 1:PARAM.nbpanels
        
        if PARAM.panels(i,end)==1 || PARAM.panels(i,end)==2 || PARAM.panels(i,end)==3
            [a,b] = drawline2(PARAM.L,0,PARAM.L,PARAM.R+PARAM.L/2*sin(PARAM.theta),PARAM.panels(i,5));
            
            eval([PARAM.varlistA{i+1} '= a']);
            eval([PARAM.varlistB{i+1} '= b']);
            
        elseif PARAM.panels(i,end)==4
            if PARAM.panels(i,3) >= PARAM.R
                b = 0.9*R;
                a = R_bubble^3/b^2;
                D = (a-b)/(a+b);
            else
                D = 0;
            end      
            
            [a,b] = draw_bubble(PARAM.L/2,0,PARAM.panels(i,5),D,PARAM.panels(i,3));
            
            eval([PARAM.varlistA{i+1} '= a']);
            eval([PARAM.varlistB{i+1} '= b']);
            
            
            
        end
        
    end

end