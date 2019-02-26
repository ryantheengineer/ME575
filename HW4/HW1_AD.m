function [F, J] = HW1_AD(x1,x2,x3,x4,x5,x6)
    
    % Other analysis variables
    h0 = 1.0;               % preload height, in
    delta0 = 0.4;           % deflection, in
    hdef = h0 - delta0;     % deflection height, in
    G = 12*10^6;            % psi
    Se = 45000;             % psi
    w = 0.18;
    Sf = 1.5;
    Q = 150000;             % psi
    delta_f = 0;            % in

    % Design variables
    d =  valder(x1,[1,0,0,0]);      % wire diameter
    D =  valder(x2,[0,1,0,0]);      % coil diameter
    n =  valder(x3,[0,0,1,0]);      % number of coils
    hf = valder(x4,[0,0,0,1]);      % unloaded height

    % Intermediate functions
    k = (G*d^4)/(8*(D^3)*n);            % spring stiffness
    delta_p = hf-h0;                    % deflection at preload
    delta_def = delta0 + (hf-h0);       % greatest working deflection
    hs = n*d;                           % solid height
    delta_s = hf - hs;                  % deflection at solid height
    F_f = k*delta_f;                    % full height force (zero)
    F_p = k*delta_p;                    % preload force
    F_def = k*delta_def;                % force at full deflection
    F_s = k*delta_s;                    % force at solid height
    K = ((4*D)-d)/(4*(D-d))+0.62*(d/D); % Wahl factor
    tau_f = ((8*F_f*D)/(pi*d^3))*K;     % stress at full height (zero)
    tau_p = ((8*F_p*D)/(pi*d^3))*K;     % stress at preload height
    tau_def = ((8*F_def*D)/(pi*d^3))*K; % stress at full deflection
    tau_s = ((8*F_s*D)/(pi*d^3))*K;     % stress at solid height
    tau_max = tau_def;                  % max stress (assumed at max deflection)
    tau_min = tau_p;                    % min stress (assumed at preload deflection)
    tau_m = (tau_max + tau_min)/2;      % mean stress
    tau_a = (tau_max - tau_min)/2;      % alternating stress
    Sy = 0.44*(Q/(d^w));                % yield stress
    dratio = D/d;                       % ratio of diameters
    dsum = D + d;                       % sum of diameters
    clash = hdef - hs;                  % clash allowance
    Sefratio = Se/Sf;                   % endurance limit to safety factor
    Syfratio = Sy/Sf;                   % yield strength to safety factor
    tau_amsum = tau_a + tau_m;          % sum of alternating and mean stresses
    
    F = [tau_a.val, tau_amsum.val, dratio.val, dsum.val, clash.val, tau_s.val];
    J = [tau_a.der, tau_amsum.der, dratio.der, dsum.der, clash.der, tau_s.der];
    
end