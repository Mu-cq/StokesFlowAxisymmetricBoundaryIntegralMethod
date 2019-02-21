%compute normal vector

function [tx,ty] = computeTangentVector(x,y,orderGeometry)

if orderGeometry==0   %straight element
   
    no = sqrt(diff(x).^2+diff(y).^2);
    nx = diff(y)./no;
    ny = -diff(x)./no;
    tx = -ny;
    ty = nx;
    
elseif orderGeometry==1   %curved element
    
    error('Not implemented')
    
end