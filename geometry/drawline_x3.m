%draw a line having (x1,y1) as starting point and (x2,y2) as ending point
%and composed by N element

function [X,Y,eta,beta]=drawline_x3(x1,y1,x2,y2,N)

%create the empty matrix
eta=zeros(1,N+1);
beta=zeros(1,N+1);
X=zeros(1,N+1);
Y=zeros(1,N+1);

%x length and y length of the line
m=x2-x1;
n=y2-y1;

%fill the horizontal vectors with the coordinates of each node and
%distribuite them with exp
for i=1:N+1
    eta(i)=(i-1)*(m/N);
    beta(i)=(i-1)*(n/N);
    if x1==x2
        X(i)=x1;
    else
        %X(i) = (eta(i)/m)^2*m+x1;
        X(i)=(exp(eta(i)*3/m)-1)*m/(exp(3)-1)+x1;
    end
    if y1==y2
        Y(i)=y1;
    else
        %Y(i) = (beta(i)/n)^1.01*n+y1;
        Y(i)=(exp(beta(i)*3/n)-1)*n/(exp(3)-1)+y1;
    end
end

% figure
% plot(X,Y,'o')
end