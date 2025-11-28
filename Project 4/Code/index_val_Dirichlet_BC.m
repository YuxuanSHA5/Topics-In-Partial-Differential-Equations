function Dbc = index_val_Dirichlet_BC(Pb,Tb,g)
% index_val_Dirichlet_BC
% 
% Dbc(:,1) = nodes
% Dbc(:,2) =  Dirichlet 

Nn = size(Pb,1);

%% 
edges = [Tb(:,[1 2]);
         Tb(:,[2 3]);
         Tb(:,[3 1])];
edges_sorted = sort(edges,2);
[uniqEdges,~,ic] = unique(edges_sorted,'rows');
counts = accumarray(ic,1);
bndEdgeIdx = find(counts==1);
bndEdges   = uniqEdges(bndEdgeIdx,:);
bndNodes   = unique(bndEdges(:));

xb = Pb(bndNodes,1);
yb = Pb(bndNodes,2);

%% Dirichlet
toly = 1e-3;

%  interface 
yI1 = 0.4;
yI2 = 0.6;

%
x1  = 0.5;
x2  = 1.0;
x3  = 1.5;

% Γ1: x≈0.5，y  (0.4,1)
g1mask  = (abs(xb - x1) < 0.03) & (yb > yI1+toly) & (yb < 1.0-toly);
g1Nodes = bndNodes(g1mask);

% Γ2: x≈1.0，y  (0,0.6)
g2mask  = (abs(xb - x2) < 0.03) & (yb > 0+toly) & (yb < yI2-toly);
g2Nodes = bndNodes(g2mask);

% Γ3: x≈1.5，y 在 (0.4,1)
g3mask  = (abs(xb - x3) < 0.03) & (yb > yI1+toly) & (yb < 1.0-toly);
g3Nodes = bndNodes(g3mask);

% ΓI: interface
gImask  = (abs(yb - yI1) < 1e-2) | (abs(yb - yI2) < 1e-2);
gINodes = bndNodes(gImask);

% ΓW: y  (0.4,0.6)
xRight  = max(xb);
gWmask  = (xb > xRight - 0.05) & (yb > yI1-toly) & (yb < yI2+toly);
gWNodes = bndNodes(gWmask);

dirNodes = [g1Nodes; g2Nodes; g3Nodes; gWNodes; gINodes];
dirVals  = [g.u1*ones(numel(g1Nodes),1);
            g.u2*ones(numel(g2Nodes),1);
            g.u3*ones(numel(g3Nodes),1);
            zeros(numel(gWNodes),1);   % ΓW = 0
            zeros(numel(gINodes),1)];  % ΓI = 0

Dbc = [dirNodes, dirVals];

end
