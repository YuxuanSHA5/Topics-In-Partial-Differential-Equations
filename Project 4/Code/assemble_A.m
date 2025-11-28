function As = assemble_A(Pb,Tb,gauss,weight,p,alpha,v)
% assemble_A  A（扩散 + 对流）
% A_ij = ∫_Omega ( alpha ∇phi_j·∇phi_i + v·∇phi_j phi_i ) dx

Nn = size(Pb,1);
Ne = size(Tb,1);
As = sparse(Nn,Nn);

nGauss = size(gauss,1);

for e = 1:Ne
    nodes  = Tb(e,:);         
    coords = Pb(nodes,:);     
    x1 = coords(1,1); y1 = coords(1,2);
    x2 = coords(2,1); y2 = coords(2,2);
    x3 = coords(3,1); y3 = coords(3,2);

    % Jacobian 
    J = [x2-x1, x3-x1;
         y2-y1, y3-y1];
    detJ = abs(det(J));
    Area = detJ/2;

   
    % hat_grad = [dphi1/dxi dphi1/deta; ...]
    hat_grad = [-1, -1;
                 1,  0;
                 0,  1];      % 3x2

    % grad_phi = J^{-T} * hat_grad
    grad_phi = (J'\hat_grad')';   % 3x2

    % 
    Ke = alpha * (grad_phi*grad_phi.') * Area;

    % 
    Ce = zeros(3,3);
    vLocal = v(nodes,:);          

    for k = 1:nGauss
        xi  = gauss(k,1);
        eta = gauss(k,2);
        w   = weight(k);

        % 
        N = [1 - xi - eta, xi, eta];   

        % 
        xg = N*coords(:,1);
        yg = N*coords(:,2); %#ok<NASGU>  

        % 
        vg = N*vLocal;  

        for j = 1:3
            vDotGrad = vg * grad_phi(j,:).';  % v·∇phi_j
            for i = 1:3
                phi_i = N(i);
                Ce(i,j) = Ce(i,j) + w * vDotGrad * phi_i * detJ;
            end
        end
    end

    Ae = Ke + Ce;
    As(nodes,nodes) = As(nodes,nodes) + Ae;
end
end
