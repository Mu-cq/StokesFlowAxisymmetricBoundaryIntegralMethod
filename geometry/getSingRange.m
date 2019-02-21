%get range

function [start,finish] = getSingRange(i,nSing)

if i==1
    start = 1;
else
    start = 1+sum(nSing(1:i-1));
end
finish = sum(nSing(1:i));