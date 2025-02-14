close all
clc
clear



% ======================================================================= %
%          Numerical Modeling of an Overrunning Clutch                %
% ======================================================================= %

% ----------------------------------------------------------------------- %
%     System description
% ----------------------------------------------------------------------- %

[ClutchPar, ShaftPar] = reading_inputs();


%...........Linear Structural MAtrices..................
NDOF = (4 + 1) * 2;         % Degrees of freedom in the system
[SMatrices.Z, SMatrices.Inertia, SMatrices.K, SMatrices.C] = Build_StructuralM(ClutchPar, ShaftPar, NDOF/2);

%% Input and output torques
Torque.Cm_max = 50;     % Maximum input torque (Nm)
Torque.tv = 0.001;      % Time for torque ramp-up (s)
ShaftPar.Cr = Torque.Cm_max;     % Resistive torque (Nm)
Fexternal = zeros(NDOF/2,1);
Fexternal(4,1) = Torque.Cm_max;

%% Simulation parameters
initial_conditions = zeros(NDOF, 1)'; % Initial state vector
v0 = 50;                   % Initial angular velocity (rad/s)
initial_conditions(8) = v0;
initial_conditions(9) = v0;


% Initialize variables
Initial.state = 'STATE1'; % Initial state of the system
Initial.Cop = 0;          % Initial torque transmitted by the clutch
Initial.C_ig = 1;         % Torque on the inner race of the clutch (BI to sprag)
Initial.C_ge = -1;        % Torque on the outer race of the clutch (sprag to BE)
Initial.NBI = 1;          % Normal force on the inner ring
Initial.TBI = 2;          % Tangential force on the inner ring

[Y] = Time_stepping(ClutchPar, ShaftPar, initial_conditions, Torque, Initial, SMatrices, Fexternal);

