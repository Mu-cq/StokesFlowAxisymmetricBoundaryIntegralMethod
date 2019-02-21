%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y]=drawline2(x1,y1,x2,y2,N)

    %create the empty matrix
    X=zeros(1,N+1);
    Y=zeros(1,N+1);

    %x length and y length of the line
    m=x2-x1;
    n=y2-y1;

    %fill the horizontal vectors with the coordinates of each node
    for i=1:N+1
        X(i)=(i-1)*(m/N)+x1;
        Y(i)=(i-1)*(n/N)+y1;
    end

end