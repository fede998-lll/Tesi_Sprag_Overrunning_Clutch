function theta_BI_BE = objective_function(N_BI, R_i, R_e, R_ge, R_gi, a)
    % objective_function - Computes the relative angular displacement (theta_BI_BE)
    %                      between the inner and outer races based on the normal force.
    %
    % Inputs:
    %   N_BI - Normal force on the inner race [N]
    %   R_i  - Radius of the inner race [m]
    %   R_e  - Radius of the outer race [m]
    %   R_ge - Radius of curvature of the sprag's outer surface [m]
    %   R_gi - Radius of curvature of the sprag's inner surface [m]
    %   a    - Distance between the centers of curvature of the sprag [m]
    %
    % Output:
    %   theta_BI_BE - Relative angular displacement between the inner and outer races [rad]

    % Parameters
    Cst = 1.57e-10;        % Hertzian contact constant
    l   = 6e-3;            % Contact length [m]
    KBI = 4500 / (15e-6);  % Stiffness of the inner race [N/m]
    KBE = 4500 / (26e-6);  % Stiffness of the outer race [N/m]
    Kg  = 4500 / (6e-6);   % Stiffness of the sprag [N/m]

    % Compute deformations
    d_BI = N_BI / KBI;                      % Elastic deformation of the inner race [m]
    d_BE = N_BI / KBE;                      % Elastic deformation of the outer race [m]
    d_H  = Cst * (N_BI^0.9 / l^0.8);        % Hertzian deformation [m]
    d_g  = N_BI / Kg;                       % Elastic deformation of the sprag [m]

    % Radii after deformation for inner and outer races
    OA  = R_i - d_BI - d_H / 2;             % Updated radius of the inner race [m]
    OA0 = R_i;                              % Initial radius of the inner race [m]
    OCi = R_i + R_gi - d_BI - d_g - d_H;    % Inner curvature radius of the sprag [m]
    OCi0 = R_i + R_gi;                      % Initial inner curvature radius [m]
    OCe = R_e - R_ge + d_BE + d_g + d_H;    % Outer curvature radius of the sprag [m]
    OCe0 = R_e - R_ge;                      % Initial outer curvature radius [m]

    % Calculate angles and projections for sprag geometry
    beta = acos((OCi^2 + OCe^2 - a^2) / (2 * OCi * OCe));   % Deformed angle between OCi and OCe [rad]
    beta0 = acos((OCi0^2 + OCe0^2 - a^2) / (2 * OCi0 * OCe0)); % Initial angle between OCi and OCe [rad]
    h = OCe * sin(beta);                   % Height of the triangle with deformation [m]
    h0 = OCe0 * sin(beta0);                % Initial height of the triangle [m]
    b = OCe * cos(beta) - OA;              % Base of the triangle after deformation [m]

    % Calculate the initial (theta_i + alpha) components
    theta_i1 = acos(h0 / OCe0);            % Angle <OCe-h
    theta_i2 = acos(h0 / a);               % Angle <hCeCi
    theta_i = theta_i1 + theta_i2;         % Initial angle theta_i (no deformation, no alpha)

    % With deformation, calculate theta_i + alpha
    theta_ia = acos((a^2 + OCe^2 - OCi^2) / (2 * a * OCe)); % Deformed angle theta_i + alpha [rad]

    % Calculate the sprag rotational angle (alpha)
    alpha = theta_ia - theta_i;

    % Compute the relative angular displacement between inner and outer races
    theta_BI_BE = -((b / (R_i - d_BI - d_H / 2)) + ...
                    ((R_ge - d_g - d_H / 2) / (R_e + d_BE + d_H / 2))) * alpha;

end

