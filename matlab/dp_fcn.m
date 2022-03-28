clc
clear all 
close all

%% symbolic
syms theta1 theta2
syms omega1 omega2
syms alpha1 alpha2
syms m1 m2 I1 I2 g
syms c1 c2 l

%% position vector --- uisng homogeneous transformation
cos1 = simplify(cos(3*pi/2 + theta1));
sin1 = simplify(sin(3*pi/2 + theta1));
R01 = [cos1,-sin1; sin1,cos1];
O01 = [0;0];            % frame1 coodinates in frame0
H01 = [R01,O01;0,0,1];

cos2 = cos(theta2);
sin2 = sin(theta2);
R12 = [cos2, -sin2; sin2, cos2];
O12 = [l;0];            % frame2 coordinates in frame1
H12 = [R12,O12;0,0,1];

G1 = H01*[c1 0 1]';
G2 = H01*H12*[c2 0 1]';
x_G1 = G1(1);           % global thigh com position x
y_G1 = G1(2);           % global thigh com position y
x_G2 = G2(1);           % global calf com position x
y_G2 = G2(2);           % global calf com position y

%% linear velocity vector --- using jacobian matrix
% easier way to find end-effector linear velocity --- instead using the whole jacobian matrix
v_G1_x = jacobian(x_G1,[theta1 theta2])*[omega1 omega2]';
v_G1_y = jacobian(y_G1,[theta1 theta2])*[omega1 omega2]';
v_G2_x = jacobian(x_G2,[theta1 theta2])*[omega1 omega2]';
v_G2_y = jacobian(y_G2,[theta1 theta2])*[omega1 omega2]';
v_G1 = [v_G1_x; v_G1_y];
v_G2 = [v_G2_x; v_G2_y];

%% angular velocity vector --- using jacobian matrix
R01_xyz = [cos1,-sin1,0; sin1,cos1,0; 0,0,1];
R12_xyz = [cos2,-sin2,0; sin2,cos2,0; 0,0,1];
R02_xyz = R01_xyz*R12_xyz;
% how to find J_w reference to (https://automaticaddison.com/the-ultimate-guide-to-jacobian-matrices-for-robotics/)
J_w = [R01_xyz*[0;0;1], R02_xyz*[0;0;1]];
w = J_w*[omega1;omega2];
disp('-> joint2 global angular velocity can be proved as omega1+omega2');
disp(['w= ',char(w)]);
disp(' ');

%% lagrangian --- find equation of motion
% I2 --- inertia respect to ground frame   (omega1+omega2) --- joint2 global angular velocity
T = 0.5*m1*(v_G1'*v_G1) + 0.5*m2*(v_G2'*v_G2) + 0.5*I1*omega1*omega1 + 0.5*I2*(omega1+omega2)*(omega1+omega2); 
V = m1*g*y_G1 + m2*g*y_G2;
L = T-V;

disp('-> copy paste energy in the code --- used for energy checking');
disp(['KE(i) = ',char(T),';']);
disp(['PE(i) = ',char(V),';']);
disp(['TE(i) = KE(i)+PE(i);']);
disp(' ');

q = [theta1 theta2];
qdot = [omega1 omega2];
qddot = [alpha1 alpha2];

for ii=1:2
    dLdqdot(ii) = diff(L,qdot(ii));
    ddt_dLdqdot(ii) = diff(dLdqdot(ii),q(1))*qdot(1) + diff(dLdqdot(ii),qdot(1))*qddot(1)+...
                      diff(dLdqdot(ii),q(2))*qdot(2) + diff(dLdqdot(ii),qdot(2))*qddot(2);
    dLdq(ii) = diff(L,q(ii));
    EOM(ii) = ddt_dLdqdot(ii) - dLdq(ii);
end

%{
    derive mass_inertia matrix
           N matrix
           G matrix
           C matrix

    how we get the mass_inertia matrix? --- notice that the difference between mass_inertia matrix and the inertia tensor
        from the format of manipulator equation of motion: M(q)qddot + c(q,qdot) + f = 0 --- M(q)qddot = F_external
        from the definition of lagrangian: ddt_dLdqdot - dLdq = F_external
        thus we have --- M(q)qddot = ddt_dLdqdot - dLdq
        according to newton's second law --- F=ma
        we can have the mass_inertia matrix M(q) = (ddt_dLdqdot - dLdq)/qddot
%}
M = jacobian(EOM,[alpha1 alpha2]);

% N is the external_force matrix without external_force since we set acceleration=0
N(1,1) = subs(EOM(1),[alpha1 alpha2],[0 0]);
N(2,1) = subs(EOM(2),[alpha1 alpha2],[0 0]);

% G is the external_force matrix only with gravity since we set acceleration = velocity = 0
G(1,1) = subs(N(1,1),[omega1 omega2],[0 0]);
G(2,1) = subs(N(2,1),[omega1 omega2],[0 0]);

% C is the external_force matrix only with coriolis force
C(1,1) = N(1,1) - G(1,1);
C(2,1) = N(2,1) - G(2,1);

disp('-> copy paste in MATLAB --- used for ode45 iteration');
disp(['M11 = ', char(M(1,1)), ';'])
disp(['M12 = ', char(M(1,2)), ';'])
disp(['M21 = ', char(M(2,1)), ';'])
disp(['M22 = ', char(M(2,2)), ';'])
disp(['C1 = ', char(C(1,1)), ';'])
disp(['C2 = ', char(C(2,1)), ';'])
disp(['G1 = ', char(G(1,1)), ';'])
disp(['G2 = ', char(G(2,1)), ';'])
disp('M = [M11 M12; M21 M22];');
disp('C = [C1; C2];');
disp('G = [G1; G2];');
disp('thetaddot = M\(-G-C);');