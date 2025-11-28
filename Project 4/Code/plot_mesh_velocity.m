%% plot_mesh_velocity.m


clear; clc; close all;

%% 1.
% mesh4.mat:  P (Nnode x 2)  
%             T (Ne x 3)     
% velocity4.mat: v (Nnode x 2)  
load('mesh4.mat');      
load('velocity4.mat'); 

node = P;        
elem = T;        
vel  = v;       

%% 2. 
figure;
triplot(elem, node(:,1), node(:,2), 'k');  
axis equal tight;
hold on;
title('Mesh with velocity field');
xlabel('x');
ylabel('y');

%% 3.  mesh 

step = 10;                            
idx  = 1:step:size(node,1);          


scale = 0.5;


quiver(node(idx,1), node(idx,2), ...  
       vel(idx,1),  vel(idx,2), ...   
       scale, 'r');                   

legend({'Mesh','Velocity arrows'}, 'Location','best');
