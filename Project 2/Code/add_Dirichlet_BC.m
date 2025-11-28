function [A, b] = add_Dirichlet_BC(A, b, Dbc)
idx = Dbc(:,1); val = Dbc(:,2);
b = b - A(:,idx)*val;
A(:,idx) = 0; A(idx,:) = 0;
for k = 1:length(idx)
    A(idx(k), idx(k)) = 1;
    b(idx(k)) = val(k);
end
end
