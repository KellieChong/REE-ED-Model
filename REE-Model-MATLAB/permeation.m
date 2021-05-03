function dCdt = permeation(t, Cions_sol, x2, D_light, D_heavy)
%permeation equations: nernst planck equation for electrochemical flux in 3d
%   Detailed explanation goes here

memSA = 0.00423;
zREE = 3;
volume = 0.0013;
I = 20; %0.03115839;
F = 96485;

% make dCdt a function of the various ion concentrations rather than
% membrane ion concentrations
% C_lightMem = K_H_light * C_H_mem ^3 *(Cions_sol(1)/Cions_sol(4) ^ 3) ^ 0.8;
% C_heavyMem = K_H_heavy * C_H_mem ^3 *(Cions_sol(2)/Cions_sol(4) ^ 3) ^ 0.8;
% C_Na_mem = 0.32 * C_H_mem *(Cions_sol(3)/Cions_sol(4));
% C_H_mem = Q_converted/3 - (C_lightMem + C_heavyMem + C_Na_mem);


commonDenom = zREE ^ 2 * D_light * x2(1) + zREE ^ 2 * D_heavy * ...
    x2(2) + 1.9^15 * x2(3) + 1.1^16 * x2(4);
    
dCdt_light = -memSA * zREE * D_light * x2(1)/commonDenom * I/(F*volume);
dCdt_heavy = -memSA * zREE * D_heavy * x2(2)/commonDenom * I/(F*volume);
dCdt_Na = -memSA * 1.9 * 10 ^ 15 * x2(3)/commonDenom * I/(F*volume);
dCdt_H = -memSA * 1.1 * 10 ^16 * x2(4)/commonDenom * I/(F*volume);

dCdt = [dCdt_light; dCdt_heavy; dCdt_Na; dCdt_H];

%syms C_light(t) C_heavy(t) C_Na(t) C_H(t)
%dCdt = [diff(C_light,t) == dCdt_light, diff(C_heavy, t) == dCdt_heavy, diff(C_Na, t) == dCdt_Na, diff(C_H, t) == dCdt_H]
   
end

