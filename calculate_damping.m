function C = calculate_damping(K1, I1, I2, zeta)
    % calculate_damping - Computes the damping coefficient C for the given system.
    %
    % Inputs:
    %   K1   - Stiffness of the shaft
    %   I1   - Inertia of the first disk
    %   I2   - Inertia of the second disk
    %   zeta - Damping ratio
    %
    % Output:
    %   C    - Damping coefficient (scalar)

    % Stiffness matrix for the system
    % Represents the stiffness coupling between the two disks
    K_matrix = [K1, -K1; 
                -K1, K1];

    % Mass matrix for the system
    % Represents the inertia properties of the two disks
    M_matrix = [I1, 0; 
                0, I2];  

    % Compute eigenvalues and eigenvectors of the system
    % The eigenvalues (lambda) correspond to the square of the natural frequencies
    % The eigenvectors (Phi) define the mode shapes
    [Phi, lambda] = eig(K_matrix, M_matrix);

    % Compute natural frequencies from eigenvalues
    omega_n = sqrt(diag(lambda)); % Natural frequencies [rad/s]

    % Compute modal damping coefficients
    % These are proportional to the natural frequencies and damping ratio
    C_mod = 2 * zeta .* sqrt(lambda); % Modal damping coefficient

    % Compute the total damping coefficient matrix
    % Use the modal matrix (Phi) to transform modal damping back to the original coordinates
    C = inv(Phi') * C_mod * inv(Phi); 

    % Extract the absolute value of the first diagonal element as the scalar damping coefficient
    % This assumes the first mode is of interest
    C = abs(C(1, 1)); 

end




