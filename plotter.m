function plotter(t, Y)
    % plotter - Plots the simulation results for angular positions and velocities.
    %
    % Inputs:
    %   t - Time vector [s]
    %   Y - Matrix of state variables:
    %       Columns:
    %       1: \theta_1   - Angular position of the input shaft
    %       2: \omega_1   - Angular velocity of the input shaft
    %       3: \theta_{BI} - Angular position of the inner race
    %       4: \omega_{BI} - Angular velocity of the inner race
    %       5: \theta_{BE} - Angular position of the outer race
    %       6: \omega_{BE} - Angular velocity of the outer race
    %       7: \theta_4   - Angular position of the output shaft
    %       8: \omega_4   - Angular velocity of the output shaft
    %       9: \theta_{sprag} - Angular position of the sprag
    %       10: \omega_{sprag} - Angular velocity of the sprag

    % Plot 1: Angular positions
    figure(1)
    hold on
    grid on
    title('Angular Positions')
    plot(t, Y(:, 1), 'r', "LineWidth", 1.5)  % \theta_1
    plot(t, Y(:, 3), 'b', "LineWidth", 1.5)  % \theta_{BI}
    plot(t, Y(:, 5), 'g', "LineWidth", 1.5)  % \theta_{BE}
    plot(t, Y(:, 7), 'k', "LineWidth", 1.5)  % \theta_4
    xlabel('Time [s]')
    ylabel('Position [rad]')
    legend('\theta_1', '\theta_{BI}', '\theta_{BE}', '\theta_4')

    % Plot 2: Angular velocities
    figure(2)
    hold on
    grid on
    title('Angular Velocities')
    plot(t, Y(:, 2), 'r', "LineWidth", 1.5)  % \omega_1
    plot(t, Y(:, 4), 'b', "LineWidth", 1.5)  % \omega_{BI}
    plot(t, Y(:, 6), 'g', "LineWidth", 1.5)  % \omega_{BE}
    plot(t, Y(:, 8), 'k', "LineWidth", 1.5)  % \omega_4
    xlabel('Time [s]')
    ylabel('Velocity [rad/s]')
    legend('\omega_1', '\omega_{BI}', '\omega_{BE}', '\omega_4')

    % Plot 3: Angular position of the sprag
    figure(3)
    hold on
    grid on
    title('Sprag Angular Position')
    plot(t, Y(:, 9), 'r', "LineWidth", 1.5)  % \theta_{sprag}
    xlabel('Time [s]')
    ylabel('Position [rad]')

    % Plot 4: Angular velocity of the sprag
    figure(4)
    hold on
    grid on
    title('Sprag Angular Velocity')
    plot(t, Y(:, 10), 'r', "LineWidth", 1.5)  % \omega_{sprag}
    xlabel('Time [s]')
    ylabel('Velocity [rad/s]')
end
