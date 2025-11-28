function phi = basis_function(p, ndx, gauss)
Ng = length(gauss);

if p == 1
    nlb = 2;
elseif p == 2
    nlb = 3;
else
    error('basis_function: p must be 1 or 2.');
end

phi = zeros(nlb, Ng);

for k = 1:Ng
    s = gauss(k);
    if p == 1
        if ndx == 0
            phi(1,k) = 1 - s;
            phi(2,k) = s;
        else
            phi(1,k) = -1;
            phi(2,k) =  1;
        end
    else
        if ndx == 0
            phi(1,k) = 2*(s-0.5).*(s-1);   % 左
            phi(2,k) = 4*s.*(1-s);         % 中
            phi(3,k) = 2*s.*(s-0.5);       % 右
        else
            phi(1,k) = 4*s - 3;
            phi(2,k) = 4 - 8*s;
            phi(3,k) = 4*s - 1;
        end
    end
end
end
