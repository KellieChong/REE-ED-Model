%constants/operating parameters
LREE = 'La'; %input('Please enter the LREE: ' , 's');
HREE = 'Nd'; %input('Please enter the HREE: ' , 's');
chelant = 'EDTA'; %input('Please enter the chelating agent: ', 's' );
C_light = 5.939311; %input('Please enter the LREE concentration in mol/L: ' ) * 1000;
C_heavy = 7.1408; %5.719555; %0.0071407773 * 1000; %input('Please enter the HREE concentration in mol/L: ' ) * 1000;
totChelantConc = 2*C_heavy; %10.39108; %0.0071407773 * 1000; %input('Please enter the chelating agent concentration in mol/L: ' ) * 1000;
pH = 5;
t = 300*60;
Q = 1.8;
n = 300;
constantH = true;

[T, Y, x1, x2] = sol(LREE, HREE, chelant, C_light, C_heavy, totChelantConc, pH, t, Q, n);
SF = Y(:, 1)./Y(:, 2);
figure()
plot(T, SF, '-o');

title('Separation Factor over time');
xlabel('Time (s)');
ylabel("Separation Factor (C LREE/C HREE)");
text(2000, SF(end)* 0.9, ['Final SF: ', num2str(SF(end))], 'FontSize', 14)

%% -------------------------------------------------------------------------------
function [T, Y, x1, x2] = sol(LREE, HREE, chelant, C_light, C_heavy, totChelantConc, pH, t, Q, n)

%for EDTA, HEDTA, and DCTA columns of the matrix, they are the K_abs values, while columns 3 and 4 are diffusion const and K_H respectively
chelant_hash = ["EDTA", "HEDTA", "DCTA", "D", "K_H"]; % the columns of the database matrix
element_hash = ["La", "Pr", "Nd", "Gd", "Y", "dissConst", "stoichiometric ratio"]; %rows of the database matrix

CA = find(chelant_hash==chelant);
LREE = find(element_hash==LREE);
HREE = find(element_hash==HREE);

%first row EDTA, 2nd row HEDTA, 3rd row DCTA
dissConst = log10([2.00, 2.67, 6.16, 10.26;...
    3.23, 5.50, 10.02, 0;...
    2.40, 3.55, 6.14, 10.26
    ]);

%map all the values to the database matrix now
database = [10^15.50, 10^13.22, 10^16.60, 4.4 * 10^13, 1.49 ; ... 
    10^16.40, 10^14.39, 10^17.01, 4.8 * 10^13, 1.40; ... 
    10^16.61, 10^14.71, 10^17.31, 5.2 * 10^13, 1.15; ... 
    10^17.37, 10^15.10, 10^18.70, 6.5 * 10^13, 0.91; ... 
    10^18.09, 10^14.49, 10^19.60, 11 * 10^13, 0.66; ... 
    0.5, 1/3, 0.5, 0, 0]; %
    
C_H = 10 ^ -pH * 1000 ;
K_abs_LREE = database(LREE, CA);
K_abs_HREE = database(HREE, CA);
dissConst = (dissConst(CA, :));
dissConst = dissConst(dissConst >= 0);

K_H_light = database(LREE, 5);
K_H_heavy = database(HREE, 5);
D_light = database(LREE, 4);
D_heavy = database(HREE, 4);
eqv = database(6,CA);
tspan = linspace (0, t, n);
Q_converted = Q*1000; %Q is now converted to mol/m3

dissEquilibrium = @(x)dissEquil(x, totChelantConc, pH, dissConst, K_abs_LREE, K_abs_HREE, C_light, C_heavy, chelant);
x0 = [2, 1, 2, 1, 0.5, 0.5];
x1 = fsolve(dissEquilibrium, x0);

%for i in (1, n) % iterate over time for changing membrane concentrations
%and append ion concentration solutions to final solution matrices

%    c0mem = [C_lightMem, C_heavyMem, C_NaMem, C_HMem] 
%    cions = [C_lightion(-1), C_heavyion(-1), C_Na(-1)]

ionEx = @(memConc)ionExchange(memConc, x1, C_H, K_H_light, K_H_heavy, eqv, Q_converted);
c0mem = [4, 4, 1, 0.1];
x2 = fsolve(ionEx, c0mem);

perm = @(t, Cions_sol)permeation(t, Cions_sol, x2, D_light, D_heavy);
y0 = [x1(1), x1(3), x1(6), C_H];
[T,Y] = ode45(perm, [0, t], y0);

Y = Y /1000; % convert back to mol/L
figure()
plot(T, Y(:, 1), '-o', T, Y(:, 2), '-o', T, Y(:, 3), '-o', T, Y(:, 4), '-o')
title('Ion Concentration over time')
xlabel("Time(s)")
ylabel("Concentration (mol/L)")
legend({'C LREE ion', 'C HREE ion', 'C Na ion', 'C H ion'}, 'location', 'best')



end

