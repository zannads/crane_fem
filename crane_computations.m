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
OMEGA = 10; % [Hz]

% compute Lk, wi


% if not ok add node in the middle or divide in third..

% once all conditions are checked 
%% create the .inp file


%% perfrom the computation required 
