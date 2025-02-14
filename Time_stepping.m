function [Y] = Time_stepping(ClutchPar,ShaftPar, initial_conditions, Torque, Initial, SMatrices, Fexternal)
global state
global C_ig
global C_ge 
    
    % Simulation time settings
    dt = 1e-5;                 % Simulation time step (s)
    t_end = 1;                 % Simulation end time (s)
    time_span = 0:dt:t_end;    % Simulation time vector
    n_plot = 0.2 / dt;         % Update plots every 0.2 seconds
    
    % Preallocate solution arrays
    NDOF = length(initial_conditions);
    num_steps = length(time_span);
    Y = zeros(num_steps, NDOF); % State variable array
    Y(1, :) = initial_conditions; % Set initial conditions
    

    C_m = @(t) Torque.Cm_max * (t / Torque.tv); % Input torque as a function of time

    %% Runge-Kutta 4th Order Simulation Loop
    for i = 1:num_steps-1
    
        % Current time and state
        t = time_span(i);
        y = Y(i, :);
        % Determine the input torque based on time
        if t < Torque.tv
            % During the ramp-up period, calculate torque using the input function
            C_m_t = C_m(t);
        else
            % After ramp-up, apply maximum torque
            C_m_t = Cm_max;
        end
        Fexternal(1,1) = C_m_t;
    
        % Compute Runge-Kutta coefficients
        dYdt1 = system_eq(t, y, ClutchPar, ShaftPar, C_m_t, Initial, SMatrices, Fexternal);
        k1 = dt * dYdt1;
        dYdt2 = system_eq(t + dt / 2, y + k1 / 2, ClutchPar, ShaftPar, C_m_t, Initial, SMatrices, Fexternal);
        k2 = dt * dYdt2;
        dYdt3 = system_eq(t + dt / 2, y + k2 / 2, ClutchPar, ShaftPar, C_m_t, Initial, SMatrices, Fexternal);
        k3 = dt * dYdt3;
        dYdt4 = system_eq(t + dt, y + k3, ClutchPar, ShaftPar, C_m_t, Initial, SMatrices, Fexternal);
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
end