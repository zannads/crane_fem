clc;
close all;
clear;

folder = cd;
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
freq_val = 10;
OMEGA_val = 2*pi*freq_val; % [Hz]

[nodes, beams] = validate_FEM(nodes, beams, Cs*OMEGA_val);
%sum( beams(:, 7) < Cs*OMEGA_val ) %should be 0
% once all conditions are checked 
%% create the .inp file
name = 'crane_ZD';
inp = '.inp';
matrices = '_mkr.mat';
frequenze = '_fre.mat';

write_inp(nodes, beams, damping, masses, [], strcat(name, inp) );

%% run the program
%estract the matrices
dmb_fem2
pause;
%% perfrom the computation required 
%load matrices
full_matrices = load( strcat(name, matrices) );

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
full_freq = load( strcat(name, frequenze) );
disp('Frequenze : ');
disp( full_freq.freq (full_freq.freq <= freq_val) );

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
vett_f=start_f:df:freq_val;
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
QF( yc_idx ,1)=1;
% xc = zeros(doc,1);
% xc(R_yO2_idx, 1) = 1;

  for k = 1:length(vett_f)
      ome = vett_f(k)*2*pi;
      A = (-ome^2*MFF+i*ome*CFF+KFF);
      % QFC = -(-ome^2*MFC+i*ome*CFC+KFC)*xc;
      % xf=A\(QF+QFC);
      % QC = (-ome^2*MCF+i*ome*CCF+KCF)*xf + (-ome^2*MCC+i*ome*CCC+KCC)*xc ;
      xf = A\QF;
      QC = (-ome^2*MCF+i*ome*CCF+KCF)*xf ;
      
      R_yO2_val = QC( R_yO2_idx );
            
      mod4(k) = abs(R_yO2_val);
      fas4(k) = angle(R_yO2_val);
  end
 
 cd(strcat(folder, '/images'));
 
 figure
 subplot 211;plot(vett_f,mod1);grid;  xlabel('freq [Hz]');  title('FRF y_A ( input Fy_A )');
 subplot 212;plot(vett_f,fas1*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 print('FRF y_A ( input Fy_A )', '-djpeg');
 figure
 subplot 211;plot(vett_f,mod2);grid; xlabel('freq [Hz]');title('FRF x_B ( input Fy_A )');
 subplot 212;plot(vett_f,fas2*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 print( 'FRF x_B ( input Fy_A )', '-djpeg');
 figure
 subplot 211;plot(vett_f,mod3);grid; xlabel('freq [Hz]');title('FRF V_{O2} ( input Fy_A )');
 subplot 212;plot(vett_f,fas3*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 print( 'FRF V_{O2} ( input Fy_A )', '-djpeg');
 figure
 subplot 211;plot(vett_f,mod4);grid;xlabel('freq [Hz]');title('FRF V_{O2} ( input Fy_C )');
 subplot 212;plot(vett_f,fas4*180/pi); grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 print( 'FRF V_{O2} ( input Fy_C )', '-djpeg');

 target_1 = max( mod3 );
 target_2 = max( mod4 );
%% QUESTION 4
xb_idx = full_matrices.idb(node_b , 1 );
xa_idx = full_matrices.idb(node_a , 1 );

Fk_L = distributed_force(0, 1, beams( 6 , end-1));
lambda6_G = [ 0, 1, 0;
             -1, 0, 0;
              0, 0, 1];
          lambda6_G = [lambda6_G, zeros(3);
                        zeros(3), lambda6_G];
Fk_G = lambda6_G'*Fk_L;

Fk_L = distributed_force(0, 1, beams( 19 , end-1));
Fk_G2 = lambda6_G'*Fk_L;

mod5 = zeros(length(vett_f), 1);
fas5 = zeros(length(vett_f), 1);
mod6 = zeros(length(vett_f), 1);
fas6 = zeros(length(vett_f), 1);
mod7 = zeros(length(vett_f), 1);
fas7 = zeros(length(vett_f), 1);

QF = zeros(dof,1);
%NODE 5
QF( full_matrices.idb(node_b , 1 ) :  full_matrices.idb(node_b , 3 ),:) = Fk_G(1:3, :);

%NODE 14
QF( full_matrices.idb(14 , 1 ) :  full_matrices.idb(14 , 3 ),:) = Fk_G(4:6, :) + Fk_G2(1:3,:);

%NODE 6
QF( full_matrices.idb(node_d , 1 ) :  full_matrices.idb(node_d , 3 ),:) = Fk_G2(4:6, :);


for k = 1:length(vett_f)
      ome = vett_f(k)*2*pi;
      A = (-ome^2*MFF+i*ome*CFF+KFF);
      xf = A\QF;
      QC = (-ome^2*MCF+i*ome*CCF+KCF)*xf ;
      
      xa_val = xf( xa_idx );
      ya_val = xf( ya_idx );
      yc_val = xf( yc_idx );
      
      mod5(k) = abs(xa_val);
      fas5(k) = angle(xa_val);
      mod6(k) = abs(ya_val);
      fas6(k) = angle(ya_val);
      mod7(k) = abs(yc_val);
      fas7(k) = angle(yc_val);
end

 figure
 subplot 211;plot(vett_f,mod5);grid;  xlabel('freq [Hz]');  title('FRF x_A ( input distribuited load )');
 subplot 212;plot(vett_f,fas5*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 print( 'FRF x_A ( input distribuited load )', '-djpeg');
 figure
 subplot 211;plot(vett_f,mod6);grid; xlabel('freq [Hz]');title('FRF y_A ( input distribuited load )');
 subplot 212;plot(vett_f,fas6*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 print( 'FRF y_A ( input distribuited load )', '-djpeg');
 figure
 subplot 211;plot(vett_f,mod7);grid; xlabel('freq [Hz]');title('FRF y_C ( input distribuited load )');
 subplot 212;plot(vett_f,fas7*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 print( 'FRF y_C ( input distribuited load )', '-djpeg');
%% QUESTION 5
% history of yA under Ma moving
vett_T = 0:0.01:T;
ma_movement = zeros(size(vett_T));
FA = zeros(size(vett_T));
ya_history = zeros(size(vett_T));
for k = 1:length(Amps)
      ome = k*2*pi/T;
      QF = zeros(dof,1);
      QF( ya_idx ,1) = MA*Amps(k)*ome^2;   %1 minus due to accelration (-ome^2) another for change of direction for acceleration 
      
      A = (-ome^2*MFF+i*ome*CFF+KFF);
      xf = A\QF;
     
      ya_val = xf( ya_idx );
      
      ma_movement = ma_movement + Amps(k)* cos( ome*vett_T + phi(k) );
      FA = FA + Amps(k)* cos( ome*vett_T + phi(k) )*MA*(ome^2); %with plus is inertia, with minus is acceleration -> remove Ma
      ya_history = ya_history + abs(ya_val)* cos( ome*vett_T + phi(k) + angle(ya_val) );
end

figure; 
%plot(vett_T, [ma_movement; A1*cos( 2*pi/T *vett_T + phi_1); A2*cos( 2*2*pi/T *vett_T + phi_2); A3*cos( 3*2*pi/T *vett_T + phi_3)] ); title('FRF yAdd (input FB)');
plot(vett_T, ma_movement ); title('y_{MA}');ylabel('m');xlabel('time [s]');
print( 'y_{MA}', '-djpeg');
figure;
plot(vett_T, FA ); title('Inertia of MA');ylabel('N');xlabel('time [s]');
print( 'Inertia of MA', '-djpeg');
figure;
plot(vett_T, ya_history); title('y_A');ylabel('m');xlabel('time [s]');
print( 'y_A', '-djpeg');
%% QUESTION 6
old_tm = sum( beams( :, 3).*beams(:, end-1) ) +MC;

%beams(4,:) = [2, 6, beams(4, 3:end) ];

% check conditions again on wi 
[nodes, beams] = validate_FEM(nodes, beams(:, 1:5), Cs*OMEGA_val);
% condition on the total mass 
new_tm = sum( beams( :, 3).*beams(:, end-1) )+MC;
 if( new_tm <= 1.05*old_tm )
     disp('First condition fullfilled');
 else
     disp('Try again');
 end

%% create the .inp file
new_name = 'crane_ZD_NEW';
write_inp(nodes, beams, damping, masses, [], strcat(new_name, inp) );
%% run the program
%estract the matrices
dmb_fem2
pause;
%% perfrom the computation required 
%load matrices
full_matrices = load( strcat(new_name, matrices) );

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
% full_freq = load( strcat(new_name, frequenze) );
% disp('Frequenze : ');
% disp( full_freq.freq (full_freq.freq <= OMEGA_val) );

mod8 = zeros(length(vett_f), 1);
fas8 = zeros(length(vett_f), 1);
mod9 = zeros(length(vett_f), 1);
fas9 = zeros(length(vett_f), 1);

QF = zeros(dof,1);
QF( ya_idx ,1)=1;

for k = 1:length(vett_f)
      ome = vett_f(k)*2*pi;
      A = (-ome^2*MFF+i*ome*CFF+KFF);
      xf = A\QF;
      QC = (-ome^2*MCF+i*ome*CCF+KCF)*xf ;
      
      %ya_val = xf( ya_idx );
      %xb_val = xf( xb_idx );
      R_yO2_val = QC( R_yO2_idx );
      
      %mod1(k) = abs(ya_val);
      %fas1(k) = angle(ya_val);
      %mod2(k) = abs(xb_val);
      %fas2(k) = angle(xb_val);
      mod8(k) = abs(R_yO2_val);
      fas8(k) = angle(R_yO2_val);
end

QF = zeros(dof,1);
QF( yc_idx ,1)=1;
% xc = zeros(doc,1);
% xc(R_yO2_idx, 1) = 1;

  for k = 1:length(vett_f)
      ome = vett_f(k)*2*pi;
      A = (-ome^2*MFF+i*ome*CFF+KFF);
      % QFC = -(-ome^2*MFC+i*ome*CFC+KFC)*xc;
      % xf=A\(QF+QFC);
      % QC = (-ome^2*MCF+i*ome*CCF+KCF)*xf + (-ome^2*MCC+i*ome*CCC+KCC)*xc ;
      xf = A\QF;
      QC = (-ome^2*MCF+i*ome*CCF+KCF)*xf ;
      
      R_yO2_val = QC( R_yO2_idx );
            
      mod9(k) = abs(R_yO2_val);
      fas9(k) = angle(R_yO2_val);
 end
 
 
 figure
 subplot 211;plot(vett_f,[mod8, mod3]);grid;  
 xlabel('freq [Hz]');  title('Difference between FRF V_{O2} ( input Fy_A )'); legend('new', 'old');
 subplot 212;plot(vett_f,[fas8, fas3]*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');
 figure
 subplot 211;plot(vett_f,[mod9, mod4]);grid; 
 xlabel('freq [Hz]');title('Difference between FRF V_{O2} ( input Fy_C )');legend('new', 'old');
 subplot 212;plot(vett_f,[fas9, fas4]*180/pi);grid; xlabel('freq [Hz]'); ylabel('degree [°]');

 if( ( max(mod8) - 0.5*target_1) <= 0 && ( max(mod9) - 0.5*target_2) <= 0 )
     disp('Request fullfilled');
 else
     disp('Try again');
 end
 
 %% FUNCTIONS
 
 function [nodes, beams] = validate_FEM(nodes,  beams, limit)
 free_node = [0, 0, 0];
 % compute Lk, wi
 idx = 1;
 n_b = size(beams, 1);
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
     if(wi <= limit)
         % half the length 4 times the freq wi... square law
         %disp(idx);
         split = ceil( sqrt(limit/wi) );
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

 end
 
 function write_inp( nodes, beams, damping, masses, springs, title)
 fileID = fopen(title,'w');
%clear evertything

 
 if( ~ isempty(nodes) )
     fprintf(fileID, '*NODES\n');
     for idx = 1:size(nodes, 1)
         fprintf(fileID, '%d \t', idx);
         fprintf(fileID,'%d %d %d \t %.1f %.1f \n', nodes(idx, :));
     end
      fprintf(fileID, '*ENDNODES\n');
 end

  if( ~ isempty(beams) )
           fprintf(fileID, '*BEAMS\n');
     for idx = 1:size(beams, 1)
         fprintf(fileID, '%d \t', idx);
         fprintf(fileID,'%d %d \t %d %.1e %.1E \n', beams(idx, 1:(end-2)) );
     end
           fprintf(fileID, '*ENDBEAMS\n');
  end
 
 if( ~ isempty(damping) )
           fprintf(fileID, '*DAMPING\n%.1f %.1e \n', damping);
 end
 
 if( ~ isempty(masses) )
           fprintf(fileID, '*MASSES\n');
     for idx = 1:size(masses, 1)
         fprintf(fileID, '%d \t', idx);
         fprintf(fileID,'%d %d %d\n', masses );
     end
           fprintf(fileID, '*ENDMASSES\n');
 end
 
 if( ~ isempty(springs) )
           fprintf(fileID, '*SPRINGS\n');
     for idx = 1:size(springs, 1)
         fprintf(fileID, '%d \t', idx);
         fprintf(fileID,'%d %d %d\n', springs );
     end
           fprintf(fileID, '*ENDSPRINGS\n');
 end
 
  fclose(fileID);
 end
 
 function Fk_L = distributed_force(px, py, Lk)
    syms e;
    fu = [(1-e/Lk); 0;0; 
          (e/Lk); 0;0];
    fv = [0; 
          2*(e/Lk)^3 - 3*(e/Lk)^2 + 1;
          Lk*( (e/Lk)^3 -2*(e/Lk)^2 + e/Lk);
          0;
          -2*(e/Lk)^3 + 3*(e/Lk)^2;
          Lk*( (e/Lk)^3 -(e/Lk)^2);];
    %only for constant px or py
    
    Fk_L = double(int(fu*px, 0, Lk) + int(fv*py, 0, Lk));
 end