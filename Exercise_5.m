%Kinematic Analysis Combo - NR Method
clear all
clc
close all
a = 0.1;
b = 0.2;
omega=1;
t = linspace(0,1,101) % looping for t 
phi = pi/6+omega*t;
u0 = [0; b + a];
J = @(u) jacobian(u, b);
eps = 1e-9;
U=[]; %Initializing empty array

for i=1:length(phi)
    
% create function handles
F = @(u) constraint(u, a, b, phi(i));
x = NR_method(F, J, u0, eps)
U=[U x];
end

theta=(U(1,:));
d=U(2,:);
figure  
plot(t, d,'r-')
title('Displacement over time');
xlabel('Time - t (s)')
ylabel('Displacement (d) (m) ')

figure 
plot(t, rad2deg(theta),'m-.')
title('Theta over time');
xlabel('Time t (s)')
ylabel('Theta (\theta) (\circ) ')

g = gradient(theta)./gradient(t) ;   

% time derivative of angle theta 
figure
plot(t,g,'r-');
title('Angular velocity over time');
xlabel('Time t (s)')
ylabel('Angular velocity $\dot{d}$ (rad/s)', 'Interpreter','latex')
g = gradient(d)./gradient(t) ;  

% time derivative of displacement 
figure 
plot(t,g,'m-.');
title('Velocity over time');
xlabel('Time-t (s)')
ylabel('Velocity $\dot{d}$ (m/s)', 'Interpreter','latex')
   

function P = constraint(u, a, b, phi)  % constraints
theta = u(1);
d = u(2);
P = [a * cos(phi) + b * cos(theta) - d;
    a * sin(phi) - b * sin(theta)];
end

function P = jacobian(u, b)   % Jacobian matrix
theta = u(1);
P = [-b * sin(theta), -1
    -b * cos(theta), 0];
end