function porous_FEM
    clear; clc; close all;

    %----------------------------------------------------
    %   -div(K grad u) = 0,  in (0,1)^2
    %   u = 0
    %   u = 1
    %   (K grad u)·n = 0
    %
    %   K = k(x,y) I
    %   k = k1
    %   k = k2
    %----------------------------------------------------

    k1 = 1.0;
    k2 = 1e-3;      %  1, 1e-1, 1e-3, 1e-5 

    h = 1/64;       
    [P, T] = generate_uniform_mesh(h);  % P: N×2,  T: Ne×3

    %  FEM 
    [A, b] = assemble_system(P, T, k1, k2);

    %  Dirichlet  (x=0 → u=0, x=1 → u=1)
    [A, b, u_D] = apply_dirichlet_bc(P, A, b);

   
    u = A \ b;

   
    u = u + u_D;
    u = max(0, min(1, u));
    
    figure;
    trisurf(T, P(:,1), P(:,2), u);
    shading interp; view(2); colorbar;
    axis equal tight;
    xlabel('x'); ylabel('y');
    title(sprintf('FEM solution u(x,y) (view2),  k_2 = %.1e', k2));

   
    nx = 200; ny = 200;
    [Xg, Yg] = meshgrid(linspace(0,1,nx), linspace(0,1,ny));

   
    F_u = scatteredInterpolant(P(:,1), P(:,2), u, 'linear', 'nearest');
    Ug  = F_u(Xg, Yg);

    figure;
    contourf(Xg, Yg, Ug, 20);   
    colorbar;
    axis equal tight;
    xlabel('x'); ylabel('y');
    title(sprintf('Contour of u(x,y),  k_2 = %.1e', k2));

  
    hx = Xg(1,2) - Xg(1,1);
    hy = Yg(2,1) - Yg(1,1);
    [Ug_y, Ug_x] = gradient(Ug, hy, hx);   % gradient(Z,y,x)

    % k(x,y)
    Kg = ones(size(Xg)) * k1;
    % x>=0.5, y>=0.5
    Kg(Xg >= 0.5 & Yg >= 0.5) = k2;
    % x<=0.5, y<=0.5
    Kg(Xg <= 0.5 & Yg <= 0.5) = k2;
   

    %  v = -K grad u
    Vx = -Kg .* Ug_x;
    Vy = -Kg .* Ug_y;

    figure;
    contourf(Xg, Yg, Ug, 20); hold on;
    colorbar;
    axis equal tight;
    xlabel('x'); ylabel('y');
    title(sprintf('Streamlines of v=-K\\nabla u (2D),  k_2 = %.1e', k2));
    streamslice(Xg, Yg, Vx, Vy);
    hold off;

    
    figure;
    trisurf(T, P(:,1), P(:,2), u);
    shading interp;
    view(3);          
    colorbar;
    axis tight;
    xlabel('x'); ylabel('y'); zlabel('u(x,y)');
    title(sprintf('3D surface of u(x,y),  k_2 = %.1e', k2));
    grid on;

  
    figure;
    contour3(Xg, Yg, Ug, 20);   
    colorbar;
    xlabel('x'); ylabel('y'); zlabel('u(x,y)');
    title(sprintf('3D contour of u(x,y),  k_2 = %.1e', k2));
    grid on;
    axis tight;

   
    figure;
    surf(Xg, Yg, Ug, 'EdgeColor','none');
    hold on;
   
    step = 5;
    Xq = Xg(1:step:end, 1:step:end);
    Yq = Yg(1:step:end, 1:step:end);
    Uq = Ug(1:step:end, 1:step:end);
    Vxq = Vx(1:step:end, 1:step:end);
    Vyq = Vy(1:step:end, 1:step:end);

   
    Vzq = zeros(size(Uq));

    quiver3(Xq, Yq, Uq, Vxq, Vyq, Vzq, 'k');

    hold off;
    xlabel('x'); ylabel('y'); zlabel('u(x,y)');
    title(sprintf('Flow field on 3D surface of u(x,y),  k_2 = %.1e', k2));
    colorbar;
    grid on;
    axis tight;
    view(3);
end

function [P, T] = generate_uniform_mesh(h)
    [X,Y] = meshgrid(0:h:1, 0:h:1);  % 规则网格
    P = [X(:), Y(:)];                % 节点坐标 N×2
    T = delaunay(P(:,1), P(:,2));    % 三角剖分 Ne×3
end


%  2.  k(x,y) 

function ke = coeff_k(xc, yc, k1, k2)
    
    if xc <= 0.5 && yc >= 0.5
        ke = k1;   
    elseif xc >= 0.5 && yc >= 0.5
        ke = k2;   
    elseif xc <= 0.5 && yc <= 0.5
        ke = k2;   
    else
        ke = k1;   
    end
end


%  3. A  b

function [A, b] = assemble_system(P, T, k1, k2)
    Np = size(P,1);      
    Ne = size(T,1);      

    A = sparse(Np, Np); 
    b = zeros(Np, 1);   

    for e = 1:Ne
     
        nodes  = T(e,:);        % 1×3
        coords = P(nodes, :);   % 3×2,  (x_i,y_i)

      
        xc = mean(coords(:,1));
        yc = mean(coords(:,2));
        ke = coeff_k(xc, yc, k1, k2);

       
        Ke = element_stiffness(coords, ke);

        A(nodes, nodes) = A(nodes, nodes) + Ke;
      
    end
end


%  4.  P1 

function Ke = element_stiffness(coords, ke)
    % coords: 3×2,  (x1,y1),(x2,y2),(x3,y3)
    x1 = coords(1,1); y1 = coords(1,2);
    x2 = coords(2,1); y2 = coords(2,2);
    x3 = coords(3,1); y3 = coords(3,2);

    %  |T|
    B = [x2-x1, x3-x1;
         y2-y1, y3-y1];
    area = abs(det(B)) / 2;

    % P1 
    % beta = [y2-y3; y3-y1; y1-y2]
    % gamma= [x3-x2; x1-x3; x2-x1]
    beta  = [y2 - y3;
             y3 - y1;
             y1 - y2];

    gamma = [x3 - x2;
             x1 - x3;
             x2 - x1];

    % grad ϕ_i = 1/(2*area) * [beta_i, gamma_i]^T
    % Ke(i,j) = ke * ∫ gradϕ_i · gradϕ_j dx
    %         = ke * area * (1/(2*area)^2) * (beta_i beta_j + gamma_i gamma_j)
    %         = ke/(4*area) * (beta*beta' + gamma*gamma')
    Ke = ke/(4*area) * (beta*beta' + gamma*gamma');
end


%  5.  Dirichlet : x=0 -> u=0,  x=1 -> u=1

function [A, b, u_D] = apply_dirichlet_bc(P, A, b)
    tol = 1e-12;
    x = P(:,1);

    left_nodes  = find(abs(x - 0) < tol);
    right_nodes = find(abs(x - 1) < tol);

    Np  = size(P,1);
    u_D = zeros(Np,1);
    u_D(right_nodes) = 1.0;    

    fixed = [left_nodes; right_nodes];

    % b = b - A*u_D
    b = b - A * u_D;

    %  u_i = u_D(i)
    for k = 1:length(fixed)
        i = fixed(k);
        A(i,:) = 0;
        A(:,i) = 0;
        A(i,i) = 1;
        b(i)   = u_D(i);
    end
end
