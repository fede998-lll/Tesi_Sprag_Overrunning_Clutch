function [ClutchPar, ShaftPar] = reading_inputs()
    %   I_1    - Inertia of the input shaft
    %   I_4    - Inertia of the output shaft
    %   I_BI   - Inertia of the inner race of the clutch
    %   I_BE   - Inertia of the outer race of the clutch
    %   I_G    - Inertia of a single sprag
    %   m_G    - Mass of a single sprag
    %   nbg    - Number of sprags
    %   Cm_max - Maximum input torque
    %   tv     - Torque ramp-up time
    %   C_m    - Input torque function handle
    %   Cr     - Resistive torque
    %   C1     - Damping coefficient for input shaft
    %   C2     - Damping coefficient for output shaft
    %   K1     - Stiffness of input shaft
    %   K2     - Stiffness of output shaft

    
    %% Input and output data for the system
    ShaftPar.I_1 = 0.6;       % Inertia of the input shaft (kg·m^2)
    ShaftPar.I_4 = 13;        % Inertia of the load shaft (kg·m^2)
    ShaftPar.K1 = 10000;      % Stiffness of the first shaft (Nm/rad)
    ShaftPar.K2 = 10000;      % Stiffness of the fourth shaft (Nm/rad)
    ShaftPar.KBI_R = 3e8;          % Radial Stiffness of the inner race [N/m]
    ShaftPar.KBE_R = 1.7308e8;     % Radial Stiffness of the outer race [N/m]
    

    
    %% Clutch parameters
    ClutchPar.R_i = 11.11e-3;   % Radius of the inner race (m)
    ClutchPar.R_e = 19.44e-3;   % Radius of the outer race (m)
    ClutchPar.R_gi = 4.02e-3;   % Radius of curvature of the sprag inner surface (m)
    ClutchPar.R_ge = 4.68e-3;   % Radius of curvature of the sprag outer surface (m)
    ClutchPar.a = 0.62e-3;      % Center-to-center distance of the sprag (m)
    ClutchPar.nbg = 12;         % Number of sprags
    ClutchPar.I_G = 1.18e-8;    % Inertia of a single sprag (kg·m^2)
    ClutchPar.m_G = 1.8e-3;     % Mass of a single sprag (kg)
    ClutchPar.mu_s = 0.09;      % Static friction coefficient 0.08<mus<0.12
    ClutchPar.mu_k = 0.06;      % Static friction coefficient mu_k = 0.5 to 0.9 mu_s 
    ClutchPar.I_BI = 0.05;     % Inertia of the inner race of the clutch
    ClutchPar.I_BE = 0.05;     % Inertia of the outer race of the clutch
    ClutchPar.length = 6e-3;         % Contact length [m]
    ClutchPar.Kg  = 7.5e8;        % Stiffness of the sprag [N/m]
    ClutchPar.Cst = 1.57e-10;     % Hertzian contact constant

    % Damping coefficients
    zeta_1 = 1/100;  % Damping ratio for the first shaft
    ShaftPar.C1 = calculate_damping(ShaftPar.K1 , ShaftPar.I_1 , ClutchPar.I_BI, zeta_1); % Damping coefficient for the input shaft (Nm·s/rad)
    zeta_2 = 1/100;  % Damping ratio for the second shaft
    ShaftPar.C2 = calculate_damping(ShaftPar.K2, ClutchPar.I_BE, ShaftPar.I_4, zeta_2); % Damping coefficient for the output shaft (Nm·s/rad)


    % Constants and stiffness values


end