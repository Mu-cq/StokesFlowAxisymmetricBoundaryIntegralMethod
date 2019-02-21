%interpolate line with splines of given coeffiecients

function out = interpolateSplines(t,a,b,c,d)

ttt = repmat(t,1,numel(a));
aaa = reshape(repmat(a,numel(t),1),1,numel(a)*numel(t));
bbb = reshape(repmat(b,numel(t),1),1,numel(b)*numel(t));
ccc = reshape(repmat(c,numel(t),1),1,numel(c)*numel(t));
ddd = reshape(repmat(d,numel(t),1),1,numel(d)*numel(t));

%splines coordinates
out = [aaa+bbb.*ttt+ccc.*ttt.^2+ddd.*ttt.^3 a(end)+b(end)+c(end)+d(end)];