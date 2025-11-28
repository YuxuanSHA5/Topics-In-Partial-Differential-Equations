function [A,b] = add_Dirichlet_BC(As,b,Dbc)
% add_Dirichlet_BC
%  As  b  Dirichlet 
%  i， u_i = g_i

A = As;

for k = 1:size(Dbc,1)
    i = Dbc(k,1);   
    g = Dbc(k,2);   

    % 
    b = b - A(:,i)*g;

    %  1
    A(i,:) = 0;
    A(:,i) = 0;
    A(i,i) = 1;

    % 
    b(i) = g;
end

end
