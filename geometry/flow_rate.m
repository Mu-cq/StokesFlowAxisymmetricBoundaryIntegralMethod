%compute flow rate given acertain velocity profile in axisymmetric
%coordinate with trapezi rule

function Q = flow_rate(r,u)

    dr = diff(r);
    int = pi*(u(1:end-1)'.*r(1:end-1)+u(2:end)'.*r(2:end)).*abs(dr);
    Q = sum(int);

end