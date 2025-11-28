function [L2err, H1err] = compute_error(u_exact, ux_exact, Pb, Tb, uh, p)
[gauss, weight] = gauss_integration(p);
phi_ref  = basis_function(p,0,gauss);
dphi_ref = basis_function(p,1,gauss);
[Ne,Nlb] = size(Tb);
Ng = length(gauss);

L2err2 = 0; H1err2 = 0;

for e = 1:Ne
    xL = Pb(Tb(e,1)); xR = Pb(Tb(e,end));
    J = xR - xL;
    for k = 1:Ng
        s = gauss(k);
        phi = phi_ref(:,k);
        dphi = dphi_ref(:,k);
        x = xL + J*s;
        uh_val = phi' * uh(Tb(e,:));
        uhx = (dphi' * uh(Tb(e,:))) / J;
        ue = u_exact(x);
        uex = ux_exact(x);
        L2err2 = L2err2 + weight(k)*(ue-uh_val)^2*abs(J);
        H1err2 = H1err2 + weight(k)*((ue-uh_val)^2 + (uex-uhx)^2)*abs(J);
    end
end

L2err = sqrt(L2err2);
H1err = sqrt(H1err2);
end
