function dYdt = system_eq(t, Y, ClutchPar, ShaftPar, C_m, Initial, SMatrices, Fexternal)

    % system_eq - Computes the derivatives of the state variables for the system.
    %
    % Inputs:
    %   t      - Current simulation time
    %   Y      - Current state vector [angles; angular velocities]

    %
    % Outputs:
    %   dYdt   - Derivative of the state vector [angular velocities; angular accelerations]



    %% Input and output data for the system
    K1 = ShaftPar.K1 ; K2 = ShaftPar.K2; 
    C1 = ShaftPar.C1; C2 = ShaftPar.C2; 
   

    % Extract angles and angular velocities from the state vector
    theta = Y(1:length(Y)/2); % Angular positions
    omega = Y(length(Y)/2+1:end); % Angular velocities

    % Calculate torques for the input and output sides
    Ce_BI = K1 * (theta(1) - theta(2)) + C1 * (omega(1) - omega(2)); % Torque on the inner race
    Ce_BE = -(K2 * (theta(3) - theta(4)) + C2 * (omega(3) - omega(4))); % Torque on the outer race

    % Calculate derivatives for the clutch subsystem
    [FNL] = ...
        clutch( Ce_BI, Ce_BE, theta, omega, C_m, ClutchPar, ShaftPar, Initial, t);

    dYdt = SMatrices.Z * Y' + [zeros(length(Y)/2,1); SMatrices.Inertia\(Fexternal + FNL)];

    % The output is a vector of angular velocities and accelerations
end
