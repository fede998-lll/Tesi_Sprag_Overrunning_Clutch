function [Z, MI, K, C] = Build_StructuralM(ClutchPar, ShaftPar, NDOF)
    MI = zeros(NDOF);
    C = zeros(NDOF);
    K  = zeros(NDOF);
    MI = MI + diag([ShaftPar.I_1, ClutchPar.I_BI, ClutchPar.I_BE,  ShaftPar.I_4, ClutchPar.I_G]);%MI + diag([ShaftPar.I_1, ClutchPar.I_BI, ClutchPar.I_BE,  ShaftPar.I_4, ClutchPar.I_G, ClutchPar.m_G]);
    %MI(5, 6) = ClutchPar.I_G;
    C(1:2,1:2) = ShaftPar.C1 *[1,-1;-1,1];
    C(3:4,3:4) = ShaftPar.C2 *[1,-1;-1,1];
    K(1:2,1:2) = ShaftPar.K1 *[1,-1;-1,1];
    K(3:4,3:4) = ShaftPar.K2 *[1,-1;-1,1];
    Z = [zeros(NDOF), eye(NDOF);
        -MI\K, -MI\C];

end