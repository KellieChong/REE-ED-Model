function y = dissEquil(x, totChelantConc, pH, dissConst, K_abs_LREE, K_abs_HREE, C_light, C_heavy, chelant)
%Solve the first system of nonlinear equations using fzero
%   Since x7 is hard coded (no unpacking available in matlab compared to
%   python), this function only solves for the case of a constant pH

x7 = 10 ^-pH * 1000 ;

if strcmp(chelant, 'HEDTA')
     alphaDenom  = prod(dissConst([1:3])) + prod(dissConst([1:2])) * x7 + ...
        prod(dissConst(1) * x7 ^ 2) + x7 ^ 3;
     
else
    alphaDenom  = prod(dissConst) + prod(dissConst([1:3]))*x7 + ...
    prod(dissConst([1:2])) * x7 ^ 2 + ...
    prod(dissConst(1) * x7 ^ 3) + x7 ^ 4;
    
end

y(1) = K_abs_LREE * x(1) * x(6) - x(2);
y(2) = K_abs_HREE * x(3) * x(6) - x(4);
y(3) = C_light - x(1) - x(2);
y(4) = C_heavy - x(3) - x(4);
y(5) = totChelantConc - x(2) - x(4) - x(5);
y(6) = x(6) - prod(dissConst([1:end]))/alphaDenom * x(5);

end

