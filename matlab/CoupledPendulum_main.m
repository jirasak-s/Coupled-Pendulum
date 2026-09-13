close all; clear; clc

%% 1. Input
% 1.1 Parameters
g = 9.80665;  % Gravitational acceleration [m/s^2]
m = 1;  % Pendulum mass [kg]
L = 1;  % Pendulum length [m]
k = 10;  % Spring constant [N/m]
d = 2;  % Distance between pivots [m]
xi = 2;  % Natural length of spring [m]

% 1.2 Initial conditions
theta10 = 0;  % Initial angle of the left pendulum [rad]
theta20 = 0;  % Initial angle of the right pendulum [rad]
thetadot10 = 0;  % Initial angular velocity of the left pendulum [rad/s]
thetadot20 = 0;  % Initial angular velocity of the right pendulum [rad/s]

% 1.3 Torque input
tau1 = 1;  % Torque at the left pendulum [N.m]
tau2 = 0;  % Torque at the right pendulum [N.m]

w = 2;  % Input frequency [rad/s]
tau = @(t) sin(w*t);  % Forcing function

%% 2. Simulation setup
% 2.1 Time domain simulation
tstart = 0; % Simulation start time [s]
tstop = 50; % Simulation end time [s]
tstep = 1e-3; % Time step [s]

tspan = tstart:tstep:tstop; % Simulation time span [s]

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

% 2.2 Frequency sweep
wstart = 2.5; % Staring frequency [rad/s]
wstop = 6; % End frequency [rad/s]
wstep = 1e-3; % Frequency step [rad/s]

wspan = wstart:wstep:wstop; % Frequency range [rad/s]

%% 3. Model initialization
% 3.1 Parameters construction
cp_param = [g m L k d xi];
F = [tau1; tau2]; % Forcing

% 3.2 Compute mass and stiffness matrices and modal properties
[M, K, Lambda, wn, Phi, U, gamma] = CoupledPendulum_ModalProp(cp_param, F);

% 3.3 Initial conditions 
theta0 = [theta10; theta20]; % Initial positions
thetadot0 = [thetadot10; thetadot20]; % Initial velocities

ICs = [theta0; thetadot0]; % ICs in physical coordinates

eta0 = Phi'*M*theta0;
etadot0 = Phi'*M*thetadot0;

ICs_modal = [eta0; etadot0]; % ICs in modal coordinates

% 3.4 State-space model (modal)
u = @(t) F*tau(t);

A = [zeros(2) eye(2); -Lambda zeros(2)];
B = [zeros(2); eye(2)];

N = @(t) gamma*tau(t);

%% 4. Simulation
% 4.1 Linear model
% 4.1.1 Time dependent
[~, eta] = ode45(@(t, x) CoupledPendulum_LIN(t, x, A, B, N), ...
                 tspan, ICs_modal, options);

Eta = [eta(:,1) eta(:,2)];
Etadot = [eta(:,3) eta(:,4)];

theta = (Phi*Eta')';
thetadot = (Phi*Etadot')';

% 4.1.2 Frequency dependent
Theta = CoupledPendulum_Freq(wspan, wn, Phi, gamma);

% 4.2 Nonlinear model
[t, theta_nl] = ode45(@(t, x) CoupledPendulum_NL(t, x, cp_param, u), ...
                  tspan, ICs, options);

%% 5. Results
% Graphic defaults
set(groot, 'defaultTextInterpreter', 'latex', ...
           'defaultTextFontName', 'Times New Roman', ...
           'defaultAxesFontName', 'Times New Roman', ...
           'defaultLineLineWidth', 1.2)

% 5.1 Mode shapes
figure(1)
subplot(2, 1, 1)
plot([0; U(:,1); 0])
title('Mode Shape: Mode 1')
xticks([]); yticks([0 1])

subplot(2, 1, 2)
plot([0; U(:,2); 0])
title('Mode Shape: Mode 2')
xticks([]); yticks([-1 0 1])

% 5.2 Time response
% 5.2.1 Linear model
figure(2)
labels = {'$\theta_1$', '$\theta_2$', ...
          '$\dot{\theta}_1$', '$\dot{\theta}_2$'};
for i =1:4
    subplot(2, 2, i)
    if i < 3
        plot(t, theta(:,i))
    else
        plot(t, thetadot(:,i - 2))
    end
    ylabel(labels{i})
    xlabel('t')
end

ax = findall(gcf, 'Type', 'axes');
linkaxes(ax, 'x')
xlim([0 t(end)])

% 5.2.2 Comparison of the linear and nonlinear models
figure(3)
subplot(2, 1, 1)
plot(t, theta(:,1), 'b')
hold on
plot(t, theta_nl(:,1), 'r')
xlabel('t')
ylabel('$\theta_1$')
legend('Linear model', 'Nonlinear model')

subplot(2, 1, 2)
plot(t, theta(:,2), 'b')
hold on
plot(t, theta_nl(:,2), 'r')
xlabel('t')
ylabel('$\theta_2$')
legend('Linear model', 'Nonlinear model')

% 5.3 Frequency response
[~, id1] = max(wspan(wspan < wn(1)));
[~, id2] = max(wspan(wspan < wn(2)));
Theta1 = Theta(:,1);
Theta2 = Theta(:,2);
if norm(gamma) ~= 0
    if gamma(1) ~= 0
        Theta1(id1 + 1) = NaN;
        Theta2(id1 + 1) = NaN;
    end

    if gamma(2) ~= 0
        Theta1(id2 + 1) = NaN;  
        Theta2(id2 + 1) = NaN;
    end

    figure(4)
    plot(wspan, Theta1, 'b')
    hold on
    plot(wspan, Theta2, 'r')
    xline(wn(1), '--')
    xline(wn(2), '--')

    ylim([-10 10]*max(F))
    legend('\Theta_1', '\Theta_2')
    xlabel('$\omega$')

    title(sprintf('Frequency response for $\\tau = [%g\\ %g]^T$', ...
                         F(1), F(2)));
end
