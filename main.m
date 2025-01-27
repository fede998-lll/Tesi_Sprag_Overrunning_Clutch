close all
clc
clear

% Define global variables
global state
global NBI
global TBI
global C_ig
global C_ge
global Cop

% Initialize global variables
state = 'STATE1'; % Initial state of the system
Cop = 0;          % Initial torque transmitted by the clutch
C_ig = 1;         % Torque on the inner race of the clutch (BI to sprag)
C_ge = -1;        % Torque on the outer race of the clutch (sprag to BE)
NBI = 1;          % Normal force on the inner ring
TBI = 2;          % Tangential force on the inner ring

%% Input and output data for the system
I_1 = 0.6;       % Inertia of the input shaft (kg·m^2)
I_4 = 13;        % Inertia of the load shaft (kg·m^2)
K1 = 10000;      % Stiffness of the first shaft (Nm/rad)
K2 = 10000;      % Stiffness of the fourth shaft (Nm/rad)
I_BI = 0.05;     % Inertia of the inner race of the clutch
I_BE = 0.05;     % Inertia of the outer race of the clutch

% Damping coefficients
zeta_1 = 1/100;  % Damping ratio for the first shaft
C1 = calculate_damping(K1, I_1, I_BI, zeta_1); % Damping coefficient for the input shaft (Nm·s/rad)
zeta_2 = 1/100;  % Damping ratio for the second shaft
C2 = calculate_damping(K2, I_BE, I_4, zeta_2); % Damping coefficient for the output shaft (Nm·s/rad)

%% Clutch parameters
R_i = 11.11e-3;   % Radius of the inner race (m)
R_e = 19.44e-3;   % Radius of the outer race (m)
R_gi = 4.02e-3;   % Radius of curvature of the sprag inner surface (m)
R_ge = 4.68e-3;   % Radius of curvature of the sprag outer surface (m)
a = 0.62e-3;      % Center-to-center distance of the sprag (m)
nbg = 12;         % Number of sprags
I_G = 1.18e-8;    % Inertia of a single sprag (kg·m^2)
m_G = 1.8e-3;     % Mass of a single sprag (kg)
mu_s = 0.09;      % Static friction coefficient

%% Input and output torques
Cm_max = 50;     % Maximum input torque (Nm)
tv = 0.001;      % Time for torque ramp-up (s)
C_m = @(t) Cm_max * (t / tv); % Input torque as a function of time
Cr = Cm_max;     % Resistive torque (Nm)

%% Simulation parameters
DOF = (4 + 1) * 2;         % Degrees of freedom in the system
initial_conditions = zeros(DOF, 1)'; % Initial state vector
v0 = 50;                   % Initial angular velocity (rad/s)
initial_conditions(6) = v0;
initial_conditions(8) = v0;

% Simulation time settings
dt = 1e-5;                 % Simulation time step (s)
t_end = 1;                 % Simulation end time (s)
time_span = 0:dt:t_end;    % Simulation time vector
n_plot = 0.2 / dt;         % Update plots every 0.2 seconds

% Preallocate solution arrays
num_steps = length(time_span);
Y = zeros(num_steps, DOF); % State variable array
Y(1, :) = initial_conditions; % Set initial conditions
Cm = zeros(length(time_span), 1); % Torque array

%% Runge-Kutta 4th Order Simulation Loop
for i = 1:num_steps-1

    % Current time and state
    t = time_span(i);
    y = Y(i, :);

    % Compute Runge-Kutta coefficients
    dYdt1 = system_eq(t, y, I_1, I_4, I_BI, I_BE, I_G, m_G, nbg, Cm_max, tv, C_m, Cr, C1, C2, K1, K2, mu_s, R_e, R_i, R_ge, R_gi, a);
    k1 = dt * dYdt1;
    dYdt2 = system_eq(t + dt / 2, y + k1 / 2, I_1, I_4, I_BI, I_BE, I_G, m_G, nbg, Cm_max, tv, C_m, Cr, C1, C2, K1, K2, mu_s, R_e, R_i, R_ge, R_gi, a);
    k2 = dt * dYdt2;
    dYdt3 = system_eq(t + dt / 2, y + k2 / 2, I_1, I_4, I_BI, I_BE, I_G, m_G, nbg, Cm_max, tv, C_m, Cr, C1, C2, K1, K2, mu_s, R_e, R_i, R_ge, R_gi, a);
    k3 = dt * dYdt3;
    dYdt4 = system_eq(t + dt, y + k3, I_1, I_4, I_BI, I_BE, I_G, m_G, nbg, Cm_max, tv, C_m, Cr, C1, C2, K1, K2, mu_s, R_e, R_i, R_ge, R_gi, a);
    k4 = dt * dYdt4;

    % Update state
    Y(i + 1, :) = (y + (k1 + 2 * k2 + 2 * k3 + k4) / 6)';

    % Print simulation data at each time step
    fprintf('%.2e | %.2e %.2e %.2e %.2e | %.2e | %.2e %.2e %.2e %.2e | %.2e | %s | %.2e %.2e\n', ...
        t, Y(i, 1), Y(i, 3), Y(i, 5), Y(i, 7), Y(i, 9), Y(i, 2), Y(i, 4), Y(i, 6), Y(i, 8), Y(i, 10), state, C_ig, C_ge);

    % Update plots every n_plot iterations
    if mod(i, n_plot) == 0 || i == num_steps - 1
        plotter(time_span(1:i+1), Y(1:i+1, :));
        drawnow; % Refresh plots immediately
    end
end

