function [Res] = Static_Equilibrium(F, ShaftPar, ClutchPar, state)
    global geom %output
    R_i = ClutchPar.R_i; R_e = ClutchPar.R_e  ; R_gi = ClutchPar.R_gi ; R_ge = ClutchPar.R_ge; a = ClutchPar.a ; 
    Cst = ClutchPar.Cst ; % Constant
    l = ClutchPar.length ;       % Contact length
    KBI_R = ShaftPar.KBI_R;      % Inner race stiffness
    KBE_R = ShaftPar.KBE_R;  % Outer race stiffness
    Kg = ClutchPar.Kg ;     % Sprag stiffness
    

    if state == 'STATE1'
        NBI = F;
    else
        TBI = F(1,1); NBI = F(2,1);
    end

    
    % Compute geometric parameters based on NBI
    geom.d_BI = NBI / KBI_R;                                    % Deformation of the inner race
    geom.d_BE = NBI / KBE_R;                                    % Deformation of the outer race
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

    if state == 'STATE1' 
        Res = sin(geom.beta) ;
    else
        TBE = -(R_i - geom.d_BI - geom.d_H / 2)*TBI/(R_e + geom.d_BE + geom.d_H / 2);
        Res = [1, cos(geom.beta), sin(geom.beta) ;
               -R_i - geom.d_BI - geom.d_H / 2, geom.b,  geom.h ]*[TBE; TBI; NBI];
    end
end

    