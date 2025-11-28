function [Pb, Tb] = FEmesh(P, T, p)
P = P(:);
Np = length(P);
Ne = size(T,1);

if p == 1
    Pb = P;
    Tb = T;
elseif p == 2
    midPts = 0.5*(P(T(:,1)) + P(T(:,2)));
    Pb = [P; midPts];
    midIdx = (Np+1:Np+Ne)';
    Tb = [T(:,1), midIdx, T(:,2)];   % 节点顺序 [左, 中, 右]
else
    error('p must be 1 or 2');
end
end
