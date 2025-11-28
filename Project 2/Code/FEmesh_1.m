function [Pb, Tb] = FEmesh(P, T, p)
% Generate FE mesh for 1D problem
% P : (Np×1) node coordinates from mesh_generation
% T : (Ne×2) element -> node connectivity (linear mesh)
% p : polynomial degree (1: linear, 2: quadratic)
%
% Pb: FE node coordinates
% Tb: FE element connectivity

P = P(:);

if p == 1
   
    Pb = P;
    Tb = T;

elseif p == 2
    
    Np = size(P, 1);   
    Ne = size(T, 1);  

  
    midPts = 0.5 * ( P(T(:,1)) + P(T(:,2)) );  

  
    Pb = [P; midPts];

  
    midIdx = (Np+1 : Np+Ne).';

   
    Tb = [T, midIdx];

else
    error('FEmesh: only p = 1 or p = 2 is supported.');
end

end
