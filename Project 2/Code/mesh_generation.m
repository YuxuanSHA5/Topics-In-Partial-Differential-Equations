function [P, T] = mesh_generation(left, right, N1)
P = linspace(left, right, N1+1)';
T = [(1:N1)', (2:N1+1)'];
end
