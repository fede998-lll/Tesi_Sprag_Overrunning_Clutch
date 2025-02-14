function [dYdt3, dYdt4, dYdt5, dYdt6, dYdt9, dYdt10] = clutch(Ce_BI, Ce_BE, theta, omega, C_m, ClutchPar, Initial, t)

    % Declare global variables to store state information
    global prev_BIddot
    global prev_BEddot
    global state
    global NBI
    global TBI
    global C_ge
    global C_ig
    global Cop
    % clutch - Computes the derivatives of the state variables for the clutch system.
    %
    % Inputs:
    %   m_G    - Mass of a single sprag
    %   I_G    - Inertia of a single sprag
    %   I_BI   - Inertia of the inner race of the clutch
    %   I_BE   - Inertia of the outer race of the clutch
    %   Ce_BI  - Elastic torque on the inner race
    %   Ce_BE  - Elastic torque on the outer race
    %   theta  - Angular positions
    %   omega  - Angular velocities
    %   C_m    - Input torque
    %   nbg    - Number of sprags
    %   R_e    - Radius of the outer race
    %   R_i    - Radius of the inner race
    %   R_ge   - Radius of curvature of the sprag outer surface
    %   R_gi   - Radius of curvature of the sprag inner surface
    %   a      - Center-to-center distance of the sprag
    %   mu_s   - Static friction coefficient
    %
    % Outputs:
    %   dYdt3  - Derivative of angular velocity of the inner race
    %   dYdt4  - Angular acceleration of the inner race
    %   dYdt5  - Derivative of angular velocity of the outer race
    %   dYdt6  - Angular acceleration of the outer race
    %   dYdt9  - Angular velocity of the sprag
    %   dYdt10 - Angular acceleration of the sprag

    %% Clutch parameters
    R_i = ClutchPar.R_i; R_e = ClutchPar.R_e  ; R_gi =ClutchPar.R_gi ; R_ge = ClutchPar.R_ge; a = ClutchPar.a ; nbg = ClutchPar.nbg; I_G = ClutchPar.I_G ; m_G = ClutchPar.m_G; 
    mu_s = ClutchPar.mu_s ; I_BI = ClutchPar.I_BI ; I_BE = ClutchPar.I_BE ;



    if t == 0
        % Initialize variables
        state = Initial.state;   Cop = Initial.Cop; C_ig = Initial.C_ig; C_ge = Initial.C_ge; NBI = Initial.NBI; TBI = Initial.TBI;          % Tangential force on the inner ring
    end


    % Transition condition to switch to STATE2
    if abs(omega(2) - omega(3)) <= 0.001*abs(omega(2))
        state = 'STATE2';
    end

    % Commented-out code for future states
    % if abs(C_ig - C_ge) < 1e-5
    %         state = 'STATE3';
    % end

    %% Freewheeling - STATE1
    if strcmp(state, 'STATE1')
        % In STATE1, the clutch is freewheeling, and a 3DOF model is used

        % Compute the normal force and geometry for the clutch
        [NBI, geom] = computing_NBI(Cop, nbg, R_e, R_i, R_ge, R_gi, a);

        % Determine the tangential force on the inner race
        if omega(2) - omega(3) < -1e-1
            TBI = +NBI * mu_s;
        else
            TBI = -NBI * mu_s;
        end

        % Dynamics matrices
        A = Ce_BE + nbg * (R_e + geom.d_BE + geom.d_H / 2);
        B = -m_G * (R_e - R_ge + geom.d_BE + geom.d_g + geom.d_H);
        C = cos(geom.beta) * TBI + sin(geom.beta) * NBI;
        D = -(R_ge - geom.d_g - geom.d_H);
        E = geom.b * TBI + geom.h * NBI;
        F = (R_ge - geom.d_g - geom.d_H / 2) / (R_e + geom.d_BE + geom.d_H / 2);

        % Coefficients for the linear system
        A1 = I_BE - A * B;
        B1 = A * B * F;
        C1 = -A * C;
        A2 = I_G - D * B;
        B2 = I_G - I_G * F + D * B * F;
        C2 = E - D * C;

        % Solve the linear system
        A = [A1, B1; A2, B2];
        b = [C1; C2];
        x = A \ b;

        % Compute derivatives for the current state
        dYdt5 = omega(3); % Angular velocity of the outer race
        dYdt6 = x(1);     % Angular acceleration of the outer race
        prev_BEddot = dYdt6;

        dYdt9 = omega(5); % Angular velocity of the sprag
        dYdt10 = x(2);    % Angular acceleration of the sprag

        dYdt3 = omega(2); % Angular velocity of the inner race
        dYdt4 = (Ce_BI + nbg * (R_i - geom.d_BI - geom.d_H / 2) * TBI) / I_BI; % Angular acceleration of the inner race
        prev_BIddot = dYdt4;
    end

    %% STATE 2 (Still not fully operational)
    if strcmp(state, 'STATE2')
        % Target angular displacement difference between inner and outer races
        theta_BI_BE_target = theta(2) - theta(3);

        % Use optimization to determine NBI based on the objective function
        options = optimset('TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off');
        NBI = fminbnd(@(N_BI) abs(objective_function(N_BI, R_i, R_e, R_ge, R_gi, a) - theta_BI_BE_target), 0, 1e5, options);

        % Compute geometry based on NBI
        Cst = 1.57e-10; % Constant
        l = 6e-3;       % Contact length
        KBI = 3e8;      % Inner race stiffness
        KBE = 1.7308e8; % Outer race stiffness
        Kg = 7.5e8;     % Sprag stiffness
        geom.d_BI = NBI / KBI;
        geom.d_BE = NBI / KBE;
        geom.d_H = Cst * ((NBI^0.9) / (l^0.8));
        geom.d_g = NBI / Kg;

        % Geometry calculations
        OA = R_i - geom.d_BI - geom.d_H / 2;
        OCi = R_i + R_gi - geom.d_BI - geom.d_g + geom.d_H;
        OCe = R_e - R_ge + geom.d_BE + geom.d_g + geom.d_H;

        % Compute beta, h, and b
        geom.beta = acos((OCi^2 + OCe^2 - a^2) / (2 * OCi * OCe));
        geom.h = OCe * sin(geom.beta);
        geom.b = OCe * cos(geom.beta) - OA;

        % Iterate until convergence
        tol = 1e-8;
        iter = 0;
        max_iter = 1000;
        converged = false;

        while ~converged && iter < max_iter
            % Compute angular accelerations
            phi_ddot = prev_BEddot - ((R_ge - geom.d_g - geom.d_H / 2) / (R_e + geom.d_BE + geom.d_H / 2)^2);
            dYdt10 = -(prev_BIddot - prev_BEddot) / (((geom.b) / (R_i - geom.d_BI - geom.d_H / 2)) + ((R_ge - geom.d_g - geom.d_H / 2) / (R_e + geom.d_BE + geom.d_H / 2)));

            % Update forces
            TBI = (I_G * (dYdt10 + phi_ddot) - m_G * (R_ge - geom.d_g - geom.d_H / 2) * (R_e - R_ge + geom.d_BE + geom.d_g + geom.d_H) * phi_ddot - NBI + ((R_ge - geom.d_g - geom.d_H / 2) * sin(geom.beta) - geom.h)) / (geom.b + (R_ge - geom.d_g - geom.d_H / 2) * cos(geom.beta));
            TBE = (-I_G * cos(geom.beta) * (dYdt10 + phi_ddot) - m_G * geom.b * (R_e - R_ge + geom.d_BE + geom.d_g + geom.d_H) * phi_ddot + NBI * (geom.b * sin(geom.beta) - geom.h * cos(geom.beta))) / (geom.b + (R_ge - geom.d_g - geom.d_H / 2) * cos(geom.beta));

            % Update angular accelerations
            dYdt4 = (Ce_BI + nbg * (R_i - geom.d_BI - geom.d_H / 2) * TBI) / I_BI;
            dYdt6 = (Ce_BE + nbg * (R_e + geom.d_BE + geom.d_H / 2) * TBE) / I_BE;

            % Check for convergence
            if abs(dYdt4 - prev_BIddot) < tol && abs(dYdt6 - prev_BEddot) < tol
                converged = true;
            else
                prev_BIddot = dYdt4;
                prev_BEddot = dYdt6;
            end
            iter = iter + 1;
        end

        % Update state derivatives
        dYdt3 = omega(2); % Angular velocity of the inner race
        dYdt5 = omega(3); % Angular velocity of the outer race
        dYdt9 = -(omega(2) - omega(3)) / (((geom.b) / (R_i - geom.d_BI - geom.d_H / 2)) + ((R_ge - geom.d_g - geom.d_H / 2) / (R_e + geom.d_BE + geom.d_H / 2))); % Sprag angular velocity
        C_ig = -TBI * (R_i - geom.d_BI - geom.d_H / 2) * nbg;
        C_ge = TBE * (R_e + geom.d_BE + geom.d_H / 2) * nbg;

        % Transition conditions
        if abs(C_ig - C_ge) < 1e-5
            Cop = C_ig;
            if abs(TBI / NBI) > mu_s
                state = 'STATE1';
            end
        end
    end

    % %% Fully engaged - STATE3
    % 
    % if strcmp(state, 'STATE3')
    % 
    %     %So che C_ig = C_ge
    % 
    %     % Valore target di theta_BI_BE
    %     theta_BI_BE_target = theta(2) - theta(3); % [rad]
    % 
    %     options = optimset('TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off');
    %     NBI = fminbnd(@(N_BI) abs(objective_function(N_BI, R_i, R_e, R_ge, R_gi, a) - theta_BI_BE_target), 0, 5e4, options);
    % 
    %     Cst = 1.57e-10;     % Costante
    %     l   = 6e-3;         % Lunghezza contatto
    %     KBI = 3e8;          % Rigidezza pista interna
    %     KBE = 1.7308e8;     % Rigidezza pista esterna
    %     Kg  = 7.5e8;        % Rigidezza sprag
    %     d_BI = NBI / KBI;
    %     d_BE = NBI / KBE;
    %     d_H  = Cst * ((NBI^0.9)/(l^0.8));
    %     d_g  = NBI / Kg;
    % 
    %     OA  = R_i - d_BI - d_H/2;
    %     OCi = R_i + R_gi - d_BI - d_g + d_H;
    %     OCe = R_e - R_ge + d_BE + d_g + d_H;
    % 
    %     % Calcolo di beta, h e b
    %     beta = acos( (OCi^2 + OCe^2 - a^2) / (2 * OCi * OCe) );
    %     h    = OCe * sin(beta);
    %     b    = OCe * cos(beta) - OA;
    % 
    %     Cop = (nbg * h * NBI) / (((R_ge - d_g - (d_H / 2)) / (R_e + d_BE + (d_H / 2))) + (b / (R_i - d_BI - (d_H / 2))));
    % 
    % 
    %     TBI = - C / ((R_i - d_BI - d_H/2)*nbg);
    %     TBE = C / ((R_e + d_BE + d_H/2)*nbg);
    % 
    %     if abs(TBI/NBI) > mu_s
    % 
    %         state = 'STATE1';
    % 
    %     end
    % 
    %     dYdt3 = omega(1);
    %     dYdt4 = (Ce_BI + nbg*(R_i - d_BI - d_H/2)*TBI)/I_BI;
    % 
    %     dYdt5 = omega(2);
    %     dYdt6 = (Ce_BE + nbg*(R_e + d_BE + d_H/2)*TBE)/I_BE;
    % 
    %     dYdt9 = 0;
    %     dYdt10 = 0;
    % 
    % end

end
