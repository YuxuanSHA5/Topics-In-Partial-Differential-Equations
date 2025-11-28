function b = assemble_b(Pb, Tb, gauss, weight, p, f)
Np = size(Pb,1);
[Ne,Nlb] = size(Tb);
Ng = length(gauss);
phi_ref = basis_function(p,0,gauss);
b = zeros(Np,1);

for e = 1:Ne
    xL = Pb(Tb(e,1)); xR = Pb(Tb(e,end));
    J = xR - xL;
    be = zeros(Nlb,1);
    for k = 1:Ng
        x = xL + J*gauss(k);
        fx = f(x);
        for i = 1:Nlb
            be(i) = be(i) + phi_ref(i,k)*fx*weight(k)*abs(J);
        end
    end
    b(Tb(e,:)) = b(Tb(e,:)) + be;
end
end
