function FNL = clutch(Ce_BI, Ce_BE, theta, omega, C_m, ClutchPar, ShaftPar, Initial, t)

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
    %   FNL

    FNL = zeros(length(theta),1);

    %% Clutch parameters
    R_i = ClutchPar.R_i; R_e = ClutchPar.R_e  ; R_gi = ClutchPar.R_gi ; R_ge = ClutchPar.R_ge; a = ClutchPar.a ; nbg = ClutchPar.nbg; I_G = ClutchPar.I_G ; m_G = ClutchPar.m_G; 
    mu_s = ClutchPar.mu_s ; I_BI = ClutchPar.I_BI ; I_BE = ClutchPar.I_BE ; mu_k = ClutchPar.mu_k;



    if t == 0
        % Initialize variables
        state = Initial.state;   Cop = Initial.Cop; C_ig = Initial.C_ig; C_ge = Initial.C_ge; NBI = Initial.NBI; TBI = Initial.TBI;          % Tangential force on the inner ring
    end


    % Transition condition to switch to STATE2
    if abs(omega(2) - omega(3)) <= 0.001*abs(omega(2))
        state = 'STATE2';
    end

    % Commented-out code for future states
    if abs(C_ig - C_ge) < 1e-5
        state = 'STATE3';
    end

    %% Freewheeling - STATE1
    if strcmp(state, 'STATE1')
        % In STATE1, the clutch is freewheeling, and a 3DOF model is used

        % Compute the normal force and geometry for the clutch
        [TBI, TBE, NBI, geom] = Contact_Forces(Cop, state, ClutchPar, ShaftPar);
        FNL(2 , 1) = nbg * (R_i - geom.d_BI - geom.d_H / 2)* TBI;
        FNL(3 , 1) = nbg * (R_e + geom.d_BE + geom.d_H / 2)* TBE;
        FNL(4 , 1) = -I_G*(TBE + TBI*cos(geom.beta)- NBI*sin(geom.beta))/m_G/(R_e - R_ge + geom.d_BE + geom.d_H + geom.d_g)-(R_ge - geom.d_g - geom.d_H / 2)*TBE + TBI * NBI*geom.h ;
    end

    %% STATE 2 (Still not fully operational)
    if strcmp(state, 'STATE2')
        % Target angular displacement difference between inner and outer races
        theta_BI_BE_target = theta(2) - theta(3);

        % Use optimization to determine NBI based on the objective function
        options = optimset('TolX', 1e-10, 'TolFun', 1e-10, 'Display', 'off');
        NBI = fminbnd(@(N_BI) abs(objective_function(N_BI, R_i, R_e, R_ge, R_gi, a) - theta_BI_BE_target), 0, 1e5, options);

        % Compute geometry based on NBI
        Cst = ClutchPar.Cst ; % Constant
        l = ClutchPar.length ;       % Contact length
        KBI_R = ShaftPar.KBI_R;      % Inner race stiffness
        KBE_R = ShaftPar.KBE_R;  % Outer race stiffness
        Kg = ClutchPar.Kg ;     % Sprag stiffness
        geom.d_BI = NBI / KBI_R;
        geom.d_BE = NBI / KBE_R;
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
        
        FNL(2 , 1) = nbg * (R_i - geom.d_BI - geom.d_H / 2)* TBI;
        FNL(3 , 1) = nbg * (R_e + geom.d_BE + geom.d_H / 2)* TBE;
        FNL(4 , 1) = -I_G*(TBE + TBI*cos(geom.beta)- NBI*sin(geom.beta))/m_G/(R_e - R_ge + geom.d_BE + geom.d_H + geom.d_g)-(R_ge - geom.d_g - geom.d_H / 2)*TBE + TBI * NBI*geom.h ;
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

    %% Fully engaged - STATE3
    
    if strcmp(state, 'STATE3')
    
        %So che C_ig = C_ge
    

        [TBI, TBE, NBI, geom] = Contact_Forces(Cop, state, ClutchPar, ShaftPar);
    
        if abs(TBI/NBI) > mu_s
    
            state = 'STATE1';
    
        end
        FNL(2 , 1) = nbg * (R_i - geom.d_BI - geom.d_H / 2)* TBI;
        FNL(3 , 1) = nbg * (R_e + geom.d_BE + geom.d_H / 2)* TBE;
        FNL(4 , 1) = 0;
   
    end

end
