function dYdt = system_eq0(t, Y, ClutchPar, ShaftPar, C_m, Initial)

    % system_eq - Computes the derivatives of the state variables for the system.
    %
    % Inputs:
    %   t      - Current simulation time
    %   Y      - Current state vector [angles; angular velocities]
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
    %
    % Outputs:
    %   dYdt   - Derivative of the state vector [angular velocities; angular accelerations]



    %% Input and output data for the system
    I_1 = ShaftPar.I_1;    I_4 = ShaftPar.I_4; K1 = ShaftPar.K1 ; K2 = ShaftPar.K2; 
    C1 = ShaftPar.C1; C2 = ShaftPar.C2; Cr = ShaftPar.Cr;
   

    % Extract angles and angular velocities from the state vector
    theta = Y(1:2:end); % Angular positions
    omega = Y(2:2:end); % Angular velocities

    % Calculate torques for the input and output sides
    Ce_BI = K1 * (theta(1) - theta(2)) + C1 * (omega(1) - omega(2)); % Torque on the inner race
    Ce_BE = -(K2 * (theta(3) - theta(4)) + C2 * (omega(3) - omega(4))); % Torque on the outer race

    % Derivatives of angular positions are equal to angular velocities
    dYdt(1) = omega(1); % Angular velocity of input shaft
    dYdt(7) = omega(4); % Angular velocity of output shaft

    % Angular accelerations for input and output shafts
    dYdt(2) = (C_m - Ce_BI) / I_1; % Input shaft angular acceleration
    dYdt(8) = (-Cr - Ce_BE) / I_4; % Output shaft angular acceleration

    % Calculate derivatives for the clutch subsystem
    [dYdt(3), dYdt(4), dYdt(5), dYdt(6), dYdt(9), dYdt(10)] = ...
        clutch( Ce_BI, Ce_BE, theta, omega, C_m, ClutchPar, Initial, t);

    % The output is a vector of angular velocities and accelerations
end
