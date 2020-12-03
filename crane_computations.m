clc;
close all;
% create the fem system with nodes on the intersections
%% data
m1 = 312; EA1 = 8.2e9; EJ1 = 1.40e9;        red = [m1, EA1, EJ1];
m2 = 200; EA2 = 5.4e9; EJ2 = 4.5e8;         green = [m2, EA2, EJ2];
m3 = 90.0; EA3 = 2.4e9; EJ3 = 2.0e8;        blue = [m3, EA3, EJ3];
MC = 2000; JC = 320;
L1 = 30; L2 = 39; L3 = 15; L4 = 24; 
H1 = 12; H2 = 15; H3 = 10;
MA = 800; 
A1 = 0.25; A2 = 0.25; A3 = 0.15; 
phi_1 = 0; phi_2 = pi; phi_3 = pi;
T = 1.2;
alpha = 0.1;
beta = 2.0e-4;

node_a = 8;
node_b = 5;
node_c = 10;
node_d = 6;
node_o1 = 1;
node_o2 = 4;

%% matrices
% types od nodes 
free_node = [0, 0, 0];
clamped_node = [1, 1, 1];
hinge_node = [1, 1, 0];

% nodes - boundary conditions codes: x,y,theta | position: x       y   
nodes = [
        hinge_node, L3, 0; 
        free_node, L3, H3;
        free_node, L3, H2+H3;
        hinge_node, L1, 0;
        free_node, L1, H3;
        free_node, L1, H2+H3;
        free_node, L1+L2-L4, H2+H3;
        free_node, L1+L2, H2+H3;
        free_node, L1, H1+H2+H3;
        free_node, 0, H2+H3;
        ];

% beams - i-th node nr.  j-th node nr.      mass [kg/m]   EA [N]  EJ [Nm^2]
% then I'll add two columns, length Lk [m] and frequency [wi]
beams = [
        1, 2, green;
        2, 3, green;
        2, 5, green;
        3, 5, green;
        4, 5, green;
        5, 6, green;
        10, 3, red;
        3, 6, red;
        6, 7, red;
        7, 8, red;
        10, 9, blue;
        3, 9, blue;
        6, 9, blue;
        7, 9, blue;
        8, 9, blue;
        ];

% alpha and beta values to define the damping matrix
damping = [alpha, beta];

% masses - node nr. where placed mass [kg]   J [Kgm^2]
masses = [ node_c, MC, JC];

%% check Cs * Ome << wi = (pi/Lk)^2 * sqrt( EJ/M )
Cs = 2;
OMEGA_val = 10; % [Hz]

% compute Lk, wi
idx = 1;
[n_b, ~] = size(beams);
beams = [beams, zeros(n_b, 2)];
while (idx <= n_b)
   % Lk
   x1 = nodes( beams(idx, 1), end-1);
   x2 = nodes( beams(idx, 2), end-1);
   y1 = nodes( beams(idx, 1), end);
   y2 = nodes( beams(idx, 2), end);
   L_x = x2 - x1;
   L_y = y2 - y1;
   Lk = sqrt( L_x^2 + L_y^2 );
 
   % wi 
   EJi = beams(idx, 5);
   Mi = beams(idx, 3);
   wi = ((pi/ Lk)^2 )*sqrt( EJi / Mi);
   
   beams(idx, end-1) = Lk;
   beams(idx, end) = wi;
   
   % if not ok add node in the middle or divide in third..
   if(wi <= Cs*OMEGA_val)
       % half the length 4 times the freq wi... square law
       disp(idx);
       split = ceil( sqrt(Cs*OMEGA_val/wi) );
       % I add nodes in the middle based on the number of pieces I had to
       % divide it ( 1 node less than beams)
       new_nodes = zeros( split -1 , 5);
       new_beams = zeros( split -1, 7);
       % now I have to fill the data of the new nodes
       for jdx = 1:split-1
           new_nodes( jdx, :) = [free_node, x1 + L_x*jdx/split, y1 + L_y*jdx/split];
       end
       [new_nodes_position, ~] = size(nodes);
       nodes = [nodes;
                new_nodes];
       %   and then change the beams
       
       %old beams needs new final node 
       old_final_node = beams(idx, 2);
       beams(idx, 2) = new_nodes_position +1 ;
       for jdx = 1:split-1
           new_beams(jdx, :) = [ new_nodes_position+ jdx, new_nodes_position+jdx+1 , beams(idx, 3:5), zeros(1, 2)];
       end
       new_beams(end, 2) = old_final_node;
       beams = [beams; 
                new_beams];
       %redo the calculations for this beam
       idx = idx-1;
       %clear new_beams new_nodes new_nodes_position old_final_node;
   end
   [n_b, ~] = size(beams);
   idx = idx +1;
end

%sum( beams(:, 7) < Cs*OMEGA_val ) %should be 0
% once all conditions are checked 
%% create the .inp file


%% perfrom the computation required 

