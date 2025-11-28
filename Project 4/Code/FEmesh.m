function [Pb,Tb] = FEmesh(P,T,p)
% FEmesh 


if p ~= 1
    error('线性单元 p = 1');
end

Pb = P;   
Tb = T;   

end
