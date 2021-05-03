function y = ionExchange(memConc, x1, C_H, K_H_light, K_H_heavy, eqv, Q_converted)
%ionExchange solves for the second system of nonlinear equations seen in
%the takahashi paper
%   Detailed explanation goes here

C_light = x1(1);
C_heavy = x1(3);
C_Na = x1(6) / eqv ;


y(1) = memConc(1) - K_H_light * memConc(4) ^3 *(C_light/C_H ^ 3) ^ 0.8;
y(2) = memConc(2) - K_H_heavy * memConc(4) ^3 *(C_heavy/C_H ^ 3) ^ 0.8;
y(3) = memConc(3) - 0.32 * memConc(4) *(C_Na/C_H);
y(4) = Q_converted - 3 * (memConc(1) + memConc(2) + memConc(3) + memConc(4));

end

