function [gauss, weight] = gauss_integration(p)
if p == 1
    gauss  = [0.211324865405187; 0.788675134594813];
    weight = [0.5; 0.5];
elseif p == 2
    gauss  = [0.112701665379258; 0.5; 0.887298334620742];
    weight = [0.277777777777778; 0.444444444444444; 0.277777777777778];
else
    error('p must be 1 or 2');
end
end
