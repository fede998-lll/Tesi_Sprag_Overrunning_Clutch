function [TBI, TBE, NBI, geom] = Contact_Forces(C, state, ClutchPar, ShaftPar, omega)
    global geom %input
    % computing_NBI - Computes the normal force (NBI) and geometric parameters for the clutch system.
    %C, nbg, R_e, R_i, R_ge, R_gi, a
    % Inputs:
    %   C    - Input torque [Nm]
    %   nbg  - Number of sprags
    %   R_e  - Radius of the outer race [m]
    %   R_i  - Radius of the inner race [m]
    %   R_ge - Radius of curvature of the sprag's outer surface [m]
    %   R_gi - Radius of curvature of the sprag's inner surface [m]
    %   a    - Distance between the centers of curvature of the sprag [m]
    %
    % Outputs:
    %   NBI  - Normal force on the inner race [N]
    %   geom - Struct containing geometric parameters:
    %          - d_BI: Elastic deformation of the inner race [m]
    %          - d_BE: Elastic deformation of the outer race [m]
    %          - d_H : Hertzian deformation [m]
    %          - d_g : Elastic deformation of the sprag [m]
    %          - beta: Angle between contact points [rad]
    %          - h   : Height of the triangle formed by the sprag [m]
    %          - b   : Base length of the triangle [m]

    R_i = ClutchPar.R_i; R_e = ClutchPar.R_e  ; R_gi = ClutchPar.R_gi ; R_ge = ClutchPar.R_ge; a = ClutchPar.a ; nbg = ClutchPar.nbg;
    Cst = ClutchPar.Cst ; % Constant
    l = ClutchPar.length ;       % Contact length
    KBI_R = ShaftPar.KBI_R;      % Inner race stiffness
    KBE_R = ShaftPar.KBE_R;  % Outer race stiffness
    Kg = ClutchPar.Kg ;     % Sprag stiffness
    if state ~= 'STATE2'
        if state == 'STATE1'
            x0 = NBI0;
        else
            x0 = [TBI0; NBI0];
        end
        fun = @(x) Static_Equilibrium(x, ShaftPar, ClutchPar, state);
        [Sol,~,~,~,~] = fsolve(fun,x0,options);
        if state == 'STATE1'
            NBI = Sol;        
            if omega(2) - omega(3) < -1e-1
                TBI = +NBI * ClutchPar.mu_k;
            else
                TBI = -NBI * ClutchPar.mu_k;
            end
            TBE = - cos(geom.beta)*TBI;
        else
            TBI = Sol(1,1); NBI = Sol(2,1);
            TBE = -(R_i - geom.d_BI - geom.d_H / 2)*TBI/(R_e + geom.d_BE + geom.d_H / 2);
        end

    end


    if state == 'STATE2'
    % Define the function to find the root of
        f = @(NBI) calcDiff(NBI, C, ClutchPar, ShaftPar);
    
        % Initial guess range for fzero
        min = 0;            % Minimum guess for NBI
        max = 7e4;          % Maximum guess for NBI
    
        % Solve for NBI using fzero
        NBI = fzero(f, [min, max]);
    
        % Compute geometric parameters based on NBI
        geom.d_BI = NBI / KBI;                                    % Deformation of the inner race
        geom.d_BE = NBI / KBE;                                    % Deformation of the outer race
        geom.d_H  = Cst * ((NBI^0.9) / (l^0.8));                  % Hertzian deformation
        geom.d_g  = NBI / Kg;                                     % Deformation of the sprag
    
        % Compute distances for geometry
        OA  = R_i - geom.d_BI - geom.d_H / 2;                     
        OCi = R_i + R_gi - geom.d_BI - geom.d_g + geom.d_H;       
        OCe = R_e - R_ge + geom.d_BE + geom.d_g + geom.d_H;       
    
        % Calculate geometric parameters beta, h, and b
        geom.beta = acos((OCi^2 + OCe^2 - a^2) / (2 * OCi * OCe)); 
        geom.h    = OCe * sin(geom.beta);                         
        geom.b    = OCe * cos(geom.beta) - OA;   
    end
end

%---------------------------------------
% Local function to compute (C - C_NBI(NBI))
%---------------------------------------
function y = calcDiff(NBI, C, nbg, R_e, R_i, R_ge, R_gi, a, Cst, l, KBI, KBE, Kg)
    % calcDiff - Computes the difference between the input torque and the
    % torque derived from the normal force (NBI).
    %
    % Inputs:
    %   NBI  - Normal force on the inner race [N]
    %   C    - Input torque [Nm]
    %   nbg  - Number of sprags
    %   R_e  - Radius of the outer race [m]
    %   R_i  - Radius of the inner race [m]
    %   R_ge - Radius of curvature of the sprag's outer surface [m]
    %   R_gi - Radius of curvature of the sprag's inner surface [m]
    %   a    - Distance between the centers of curvature of the sprag [m]
    %   Cst, l, KBI, KBE, Kg - Constants and stiffnesses
    %
    % Output:
    %   y    - Difference between the actual and calculated torque

    % Compute deformations
    d_BI = NBI / KBI;                                       % Deformation of the inner race
    d_BE = NBI / KBE;                                       % Deformation of the outer race
    d_H  = Cst * ((NBI^0.9) / (l^0.8));                     % Hertzian deformation
    d_g  = NBI / Kg;                                        % Deformation of the sprag

    % Compute distances for geometry
    OA  = R_i - d_BI - d_H / 2;                             
    OCi = R_i + R_gi - d_BI - d_g + d_H;                    
    OCe = R_e - R_ge + d_BE + d_g + d_H;                    

    % Calculate beta, h, and b
    beta = acos((OCi^2 + OCe^2 - a^2) / (2 * OCi * OCe));    
    h    = OCe * sin(beta);                                 
    b    = OCe * cos(beta) - OA;                            

    % Calculate the torque derived from NBI
    C_NBI = (nbg * h * NBI) / (((R_ge - d_g - d_H / 2) / (R_e + d_BE + d_H / 2)) + ...
                                (b / (R_i - d_BI - d_H / 2)));

    % Compute the difference between the input torque and calculated torque
    y = C - C_NBI;
end
