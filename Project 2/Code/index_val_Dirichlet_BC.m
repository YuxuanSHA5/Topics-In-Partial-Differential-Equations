function Dbc = index_val_Dirichlet_BC(left, right, Pb, g)
tol = 1e-12;
idx_left  = find(abs(Pb-left) < tol);
idx_right = find(abs(Pb-right) < tol);
index = [idx_left; idx_right];
val = g(Pb(index));
Dbc = [index, val];
end
