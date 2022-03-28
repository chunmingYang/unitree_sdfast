clc
clear all
close all
%{
    --- ode45 iteration routine
    --- matlab animation routine
%}

%% parameters initialization
m1 = 1.013; m2 = 0.226; g = 9.8; 
c1 = 0.1; c2 = 0.1; l = 0.2;
I1 = m1*c1^2; I2 = m2*(l+c2)^2; % I1 I2 both inertia respect to ground frame
fps = 30;

parms.m1 = m1; parms.m2 = m2; parms.g = g;
parms.c1 = c1; parms.c2 = c2; parms.l = l;
parms.I1 = I1; parms.I2 = I2;

%% iteration main loop
t = linspace(0,10,1001); % time values
z0 = [90*pi/180 0 0 0];  % initialization --- theta1 omega1 theta2 omega2
options = odeset('Abstol',1e-10,'Reltol',1e-10);
[t, z] = ode45(@rhs,t,z0,options,parms); % integrate using ode45

%% result checking

% energy checking
for i = 1:length(t)
    theta1 = z(i,1); omega1 = z(i,2); theta2 = z(i,3); omega2 = z(i,4);
    KE(i) = (m1*(c1*omega1*cos(conj(theta1))*conj(c1)*conj(omega1)*cos(theta1) + c1*omega1*sin(conj(theta1))*conj(c1)*conj(omega1)*sin(theta1)))/2 + (I2*(omega1 + omega2)^2)/2 + (m2*((omega1*(cos(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))) + c2*omega2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))))*(conj(omega1)*(l*cos(theta1) + conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + conj(c2)*conj(omega2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + (omega1*(sin(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))) + c2*omega2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))))*(conj(omega1)*(l*sin(theta1) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))) + conj(c2)*conj(omega2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)))))/2 + (I1*omega1^2)/2;
    PE(i) = - g*m2*(l*cos(theta1) + conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) - g*m1*conj(c1)*cos(theta1);
    TE(i) = KE(i)+PE(i);
end

% data from sdfast
data_sd = dlmread('data_sd.csv', ',',1,0);
t_sd = data_sd(:,1);
pe_sd = data_sd(:,2);
ke_sd = data_sd(:,3);
te_sd = data_sd(:,4);
q1_sd = data_sd(:,5);
q2_sd = data_sd(:,6);

figure(1)
z_sd = [q1_sd,q1_sd,q2_sd,q2_sd];
% animate(t,z,parms,fps);  % animation_matlab
animate(t_sd,z_sd,parms,fps); % animation_sdfast

figure(2)
subplot(2,1,1);
plot(t, z(:,1),'LineWidth',3); hold on;
plot(t_sd, q1_sd,'-o')
legend('matlab','sdfast');
xlabel('time');
ylabel('theta1');
subplot(2,1,2);
plot(t, z(:,3),'LineWidth',3); hold on;
plot(t_sd, q2_sd, '-o');
legend('matlab','sdfast');
xlabel('time');
ylabel('theta2');

figure(3)
subplot(3,1,1);
plot(t,KE,'LineWidth',3); hold on;
plot(t_sd, ke_sd,'-o'); hold on;
legend('matlab','sdfast');
xlabel('time');
ylabel('KE')
subplot(3,1,2);
plot(t,PE,'LineWidth',3); hold on;
plot(t_sd,pe_sd,'-o'); hold on;
legend('matlab','sdfast');
xlabel('time');
ylabel('PE')
subplot(3,1,3)
plot(t,TE,'LineWidth',3); hold on;
plot(t_sd, te_sd); hold on;
ylim([-0.01 0.01])
legend('matlab','sdfast');
xlabel('time');
ylabel('TE')

%% iteration routine
function zdot = rhs(t,z,parms)
    m1 = parms.m1; m2 = parms.m2; g = parms.g;
    c1 = parms.c1; c2 = parms.c2; l = parms.l;
    I1 = parms.I1; I2 = parms.I2;
    theta1 = z(1);
    omega1 = z(2);
    theta2 = z(3);
    omega2 = z(4);
    tau_1 = 0;
    tau_2 = 0;
    tau = [tau_1; tau_2]; % external_force --- control torque
    
    M11 = I1 + I2 + (m2*(2*(cos(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))))*(l*cos(theta1) + conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + 2*(sin(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))))*(l*sin(theta1) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)))))/2 + (m1*(2*c1*cos(conj(theta1))*conj(c1)*cos(theta1) + 2*c1*sin(conj(theta1))*conj(c1)*sin(theta1)))/2;
    M12 = I2 + (m2*(conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))*(cos(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))*(sin(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))*(l*cos(theta1) + conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))*(l*sin(theta1) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)))))/2;
    M21 = I2 + (m2*(conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))*(cos(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))*(sin(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))*(l*cos(theta1) + conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))*(l*sin(theta1) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)))))/2;
    M22 = I2 + (m2*(2*c2*conj(c2)*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + 2*c2*conj(c2)*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))))/2;
    C1 = -(m2*omega2*((c2*omega1*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))) + c2*omega2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))))*(l*cos(theta1) + conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) - (c2*omega1*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))) + c2*omega2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))))*(l*sin(theta1) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))) + (conj(c2)*conj(omega1)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + conj(c2)*conj(omega2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)))*(cos(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))) - (conj(c2)*conj(omega1)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + conj(c2)*conj(omega2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)))*(sin(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))) + conj(c2)*(omega1*(cos(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))) + c2*omega2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))))*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) - conj(c2)*(omega1*(sin(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))) + c2*omega2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))))*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))*(conj(omega1)*(l*cos(theta1) + conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + conj(c2)*conj(omega2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) - c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))*(conj(omega1)*(l*sin(theta1) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))) + conj(c2)*conj(omega2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)))))/2;
    C2 = (m2*((omega1*(cos(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))) + c2*omega2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))))*(conj(c2)*conj(omega1)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + conj(c2)*conj(omega2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))) - (omega1*(sin(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))) + c2*omega2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))))*(conj(c2)*conj(omega1)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + conj(c2)*conj(omega2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + (c2*omega1*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))) + c2*omega2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))))*(conj(omega1)*(l*cos(theta1) + conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + conj(c2)*conj(omega2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) - (c2*omega1*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))) + c2*omega2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))))*(conj(omega1)*(l*sin(theta1) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))) + conj(c2)*conj(omega2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)))))/2 - (m2*omega2*(conj(c2)*(omega1*(cos(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))) + c2*omega2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))))*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) - conj(c2)*(omega1*(sin(conj(theta1))*conj(l) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))) + c2*omega2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))))*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))*(conj(c2)*conj(omega1)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + conj(c2)*conj(omega2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))) - c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))*(conj(c2)*conj(omega1)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + conj(c2)*conj(omega2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + c2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1)))*(conj(omega1)*(l*cos(theta1) + conj(c2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + conj(c2)*conj(omega2)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2))) + conj(c2)*(c2*omega1*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))) + c2*omega2*(cos(conj(theta1))*sin(conj(theta2)) + cos(conj(theta2))*sin(conj(theta1))))*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) - conj(c2)*(c2*omega1*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))) + c2*omega2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2))))*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) - c2*(cos(conj(theta1))*cos(conj(theta2)) - sin(conj(theta1))*sin(conj(theta2)))*(conj(omega1)*(l*sin(theta1) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))) + conj(c2)*conj(omega2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)))))/2;
    G1 = g*m2*(l*sin(theta1) + conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1))) + g*m1*conj(c1)*sin(theta1);
    G2 = g*m2*conj(c2)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1));
    M = [M11 M12; M21 M22];
    C = [C1; C2];
    G = [G1; G2];
    thetaddot = M\(tau-G-C);
    zdot = [omega1 thetaddot(1) omega2 thetaddot(2)]';
end

%% animation function
function animate(t_all, z_all, parms, fps)
    z_all_plot = [z_all(:,1), z_all(:,3)];         % animation only needs theta1 theta2
    nn = size(z_all_plot, 2);                      % number of features
    total_frames = round(t_all(end)*fps);          % fps sample step_3
    t = linspace(0,t_all(end),total_frames);       % fps sample step_2
    z = zeros(total_frames,nn);                    % fps sample step_3 
    for i=1:nn
        z(:,i) = interp1(t_all,z_all_plot(:,i),t); % matlab table look up function
    end
    l = parms.l;
    c1 = parms.c1;
    c2 = parms.c2;
    ll = 2.2*l;
    mm = size(z,1);                                % number of observations
    for i=1:mm
        theta1 = z(i,1);
        theta2 = z(i,2);
        O = [0 0];
        G1 = [c1*sin(theta1) -c1*cos(theta1)];
        P = [l*sin(theta1) -l*cos(theta1)];
        G2 = P + c2*[sin(theta1+theta2) -cos(theta1+theta2)];
        Q = P + l*[sin(theta1+theta2) -cos(theta1+theta2)];
        h1 = plot(G1(1),G1(2),'ko','MarkerFaceColor','k','Markersize',10); hold on;
        h2 = plot(G2(1),G2(2),'ko','MarkerFaceColor','k','Markersize',10);
        h3 = line([O(1) P(1)],[O(2) P(2)],'Color','red','Linewidth',2);
        h4 = line([P(1) Q(1)],[P(2) Q(2)],'Color','red','Linewidth',2);
        axis('equal')
        axis([-ll ll -ll ll]);
        if (i==1)
            pause(1)
        end
        pause(0.01);
        if (i~=mm)                                 % delete all but not the last entry
            delete(h1);
            delete(h2);
            delete(h3);
            delete(h4);
        end
   end
end





