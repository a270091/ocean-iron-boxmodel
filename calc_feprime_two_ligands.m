% test script to check how the calculation of the free iron concentration
% from total dFe and two ligands works

FeT = 1;
L1T = 1;
L2T = 1;
K1  = 100;
K2  = 1000;

% coefficients of third-order polynomial
a3 = -K1*K2;
a2 = ( K1*K2*(FeT - L1T - L2T) - K1 - K2 );
a1 = (K1 + K2)*FeT - 1 - K1*L1T - K2*L2T;
a0 = FeT;

% solve cubic equation
test = sort(roots([a3 a2 a1 a0]),'ComparisonMethod','real');

% calculate speciation
Feprime = test(end);
FeL1 = K1 * L1T * Feprime / (1 + K1*Feprime);
FeL2 = K2 * L2T * Feprime / (1 + K2*Feprime);
L1 = L1T - FeL1;
L2 = L2T - FeL2;

% print speciation and mass balances:
fprintf('Fe'':  %8.4f\n',Feprime)
fprintf('FeL1: %8.4f\n',FeL1)
fprintf('FeL2: %8.4f\n',FeL2)

% check whether speciation fits with the original equations
R1 = FeT - Feprime - FeL1 - FeL2;
R2 = K1 * L1 * Feprime - FeL1;
R3 = K2 * L2 * Feprime - FeL2;

