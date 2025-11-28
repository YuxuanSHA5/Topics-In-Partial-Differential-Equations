function [A, b] = add_Dirichlet_BC(A, b, Dbc)

% Impose Dirichlet boundary conditions on system A u = b

% Dbc: two-column matrix [node_index, boundary_value]



index = Dbc(:,1);    
val   = Dbc(:,2);    
Ni    = length(index);


b = b - A(:, index) * val;


A(:, index) = 0;
A(index, :) = 0;


b(index) = val;
for k = 1:Ni
    A(index(k), index(k)) = 1;
end

end
