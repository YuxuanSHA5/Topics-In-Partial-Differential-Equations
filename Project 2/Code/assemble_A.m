function A = assemble_A(Pb, Tb, gauss, weight, p, mu)
Np = size(Pb,1);
[Ne,Nlb] = size(Tb);
Ng = length(gauss);
dphi_ref = basis_function(p,1,gauss);
A = sparse(Np,Np);

for e = 1:Ne
    xL = Pb(Tb(e,1)); xR = Pb(Tb(e,end));
    J = xR - xL;
    Ae = zeros(Nlb,Nlb);
    for i = 1:Nlb
        for j = 1:Nlb
            for k = 1:Ng
                Ae(i,j) = Ae(i,j) + mu*(dphi_ref(i,k)/J)*(dphi_ref(j,k)/J)*weight(k)*abs(J);
            end
        end
    end
    A(Tb(e,:), Tb(e,:)) = A(Tb(e,:), Tb(e,:)) + Ae;
end
end
