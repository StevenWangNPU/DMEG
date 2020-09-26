% Given the distance to construct a sparse graph
function  S =  construct_S6(dis, idx, h, lambda, ni)
S = zeros(ni,ni);
a = zeros(ni,ni);
for j = 1:ni
    a(j) = sum(exp(-dis(j,2:h+1)./lambda));
    for k = 1:h
        S(j, idx(j,k+1)) = exp(-dis(j,k+1)./lambda)/a(j);
    end
end
