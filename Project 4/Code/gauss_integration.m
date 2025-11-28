function [gauss,weight] = gauss_integration(order)
% gauss_integration 
% (0,0), (1,0), (0,1)

% 
gauss = [1/6, 1/6;
         2/3, 1/6;
         1/6, 2/3];
weight = (1/6)*ones(3,1);  

end
