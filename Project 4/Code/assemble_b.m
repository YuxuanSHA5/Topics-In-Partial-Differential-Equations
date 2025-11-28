function b = assemble_b(Pb,Tb,gauss,weight,p,f)
% assemble_b  b
% b_i = ∫_Omega f(x,y) phi_i dx


Nn = size(Pb,1);
Ne = size(Tb,1);
b  = zeros(Nn,1);

nGauss = size(gauss,1);

for e = 1:Ne
    nodes  = Tb(e,:);
    coords = Pb(nodes,:);
    x1 = coords(1,1); y1 = coords(1,2);
    x2 = coords(2,1); y2 = coords(2,2);
    x3 = coords(3,1); y3 = coords(3,2);

    J = [x2-x1, x3-x1;
         y2-y1, y3-y1];
    detJ = abs(det(J));

    be = zeros(3,1);

    for k = 1:nGauss
        xi  = gauss(k,1);
        eta = gauss(k,2);
        w   = weight(k);

        N = [1 - xi - eta, xi, eta];    % 形函数
        xg = N*coords(:,1);
        yg = N*coords(:,2);

        fk = f(xg,yg);

        for i = 1:3
            be(i) = be(i) + w * fk * N(i) * detJ;
        end
    end

    b(nodes) = b(nodes) + be;
end
end
