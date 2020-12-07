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
Amps = [A1, A2, A3];
phi_1 = 0; phi_2 = pi; phi_3 = pi;
phi = [phi_1, phi_2, phi_3];
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
       %disp(idx);
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

%% run the program
%estract the matrices

%% perfrom the computation required 
%load matrices
full_matrices = load('crane_ZD_mkr.mat');

dof = sum (sum( nodes(:, 1:3) == [0, 0, 0]) );
doc = (3*length(nodes) - dof);
MFF = full_matrices.M( 1:dof, 1:dof);
CFF = full_matrices.R( 1:dof, 1:dof);
KFF = full_matrices.K( 1:dof, 1:dof);

MFC = full_matrices.M( 1:dof, dof+1:end);
CFC = full_matrices.R( 1:dof, dof+1:end);
KFC = full_matrices.K( 1:dof, dof+1:end);

MCF = full_matrices.M( dof+1:end, 1:dof);
CCF = full_matrices.R( dof+1:end, 1:dof);
KCF = full_matrices.K( dof+1:end, 1:dof);

MCC = full_matrices.M( dof+1:end, dof+1:end);
CCC = full_matrices.R( dof+1:end, dof+1:end);
KCC = full_matrices.K( dof+1:end, dof+1:end);
% modes and frequencies
% [modes, eigenvalues] = eig(MFF\KFF);
% freq = sqrt(diag(eigenvalues))/2/pi;
full_freq = load('crane_ZD_fre.mat');
disp('Frequenze : ');
disp( full_freq.freq (full_freq.freq <= OMEGA_val) );

%% QUESTION 3

%3.a F yA -> ya
%3.b F yA -> xb
%3.c F yA -> RyO2
%3.d F yc -> RyO2
% to resort the index of yA in the big matrices
% Ya = full_matrices.idb(node_a , 2 ) % row is the node number, column is
% the x(1), y(2) or theta(3)
ya_idx = full_matrices.idb(node_a , 2 );
xb_idx = full_matrices.idb(node_b , 1 );
yc_idx = full_matrices.idb(node_c , 2 );
R_yO2_idx = full_matrices.idb(node_o2 , 2 ) - dof ;

start_f = 0;
df = 0.001;
i=sqrt(-1);
vett_f=start_f:df:OMEGA_val;
mod1 = zeros(length(vett_f), 1);
fas1 = zeros(length(vett_f), 1);
mod2 = zeros(length(vett_f), 1);
fas2 = zeros(length(vett_f), 1);
mod3 = zeros(length(vett_f), 1);
fas3 = zeros(length(vett_f), 1);
mod4 = zeros(length(vett_f), 1);
fas4 = zeros(length(vett_f), 1);

QF = zeros(dof,1);
QF( ya_idx ,1)=1;

for k = 1:length(vett_f)
      ome = vett_f(k)*2*pi;
      A = (-ome^2*MFF+i*ome*CFF+KFF);
      xf = A\QF;
      QC = (-ome^2*MCF+i*ome*CCF+KCF)*xf ;
      
      ya_val = xf( ya_idx );
      xb_val = xf( xb_idx );
      R_yO2_val = QC( R_yO2_idx );
      
      mod1(k) = abs(ya_val);
      fas1(k) = angle(ya_val);
      mod2(k) = abs(xb_val);
      fas2(k) = angle(xb_val);
      mod3(k) = abs(R_yO2_val);
      fas3(k) = angle(R_yO2_val);
end

QF = zeros(dof,1);
xc = zeros(doc,1);
xc(R_yO2_idx, 1) = 1;

  for k = 1:length(vett_f)
      ome = vett_f(k)*2*pi;
      A = (-ome^2*MFF+i*ome*CFF+KFF);
      QFC = -(-ome^2*MFC+i*ome*CFC+KFC)*xc;
      xf=A\(QF+QFC);
      QC = (-ome^2*MCF+i*ome*CCF+KCF)*xf + (-ome^2*MCC+i*ome*CCC+KCC)*xc ;
      
      R_yO2_val = QC( R_yO2_idx );
            
      mod4(k) = abs(R_yO2_val);
      fas4(k) = angle(R_yO2_val);
 end
 
 figure
 subplot 211;plot(vett_f,mod1);grid;  xlabel('freq [Hz]');  title('FRF yA (input FA)');
 subplot 212;plot(vett_f,fas1*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 figure
 subplot 211;plot(vett_f,mod2);grid; xlabel('freq [Hz]');title('FRF yB (input FA)');
 subplot 212;plot(vett_f,fas2*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 figure
 subplot 211;plot(vett_f,mod3);grid; xlabel('freq [Hz]');title('FRF yAdd (input FB)');
 subplot 212;plot(vett_f,fas3*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 figure
 subplot 211;plot(vett_f,mod4);grid;xlabel('freq [Hz]');title('FRF yBdd (input FB)');
 subplot 212;plot(vett_f,fas4*180/pi); grid; xlabel('freq [Hz]'); ylabel('degree [°]');


%% QUESTION 4


%% QUESTION 5
% history of yA under Ma moving
vett_T = 0:0.01:T;
ma_movement = zeros(size(vett_T));
FA = zeros(size(vett_T));
ya_history = zeros(size(vett_T));
for k = 1:3
      ome = k*2*pi/T;
      QQF = zeros(dof,1);
      QF( ya_idx ,1) = MA*Amps(k)*ome^2;   %1 minus due to accelration (-ome^2) another for change of direction for acceleration 
      
      A = (-ome^2*MFF+i*ome*CFF+KFF);
      xf = A\QF;
     
      ya_val = xf( ya_idx );
      
      ma_movement = ma_movement + Amps(k)* cos( ome*vett_T + phi(k) );
      FA = FA + Amps(k)* cos( ome*vett_T + phi(k) )*MA*(ome^2);
      ya_history = ya_history + abs(ya_val)* cos( ome*vett_T + phi(k) + angle(ya_val) );
end

figure; 
%plot(vett_T, [ma_movement; A1*cos( 2*pi/T *vett_T + phi_1); A2*cos( 2*2*pi/T *vett_T + phi_2); A3*cos( 3*2*pi/T *vett_T + phi_3)] ); title('FRF yAdd (input FB)');
plot(vett_T, ma_movement ); title('FRF yAdd (input FB)');
figure;
plot(vett_T, FA ); title('FA');
figure;
plot(vett_T, ya_history);


%% QUESTION 6
