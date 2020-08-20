% test script to check how the calculation of the free iron concentration
% from total dFe and three ligands works

FeT = 1;
L1T = 1;
L2T = 1;
L3T = 1;
K1  = 100;
K2  = 1000;
K3  = 10;

%-----------
% first guess from solving with two strong ligands only 
%-----------

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
FeL3 = 0;
L1 = L1T - FeL1;
L2 = L2T - FeL2;
L3 = L3T - FeL3;

sol0 = [Feprime;FeL1;FeL2;FeL3];

%-----------
% solution for three ligands is done by Newton-Raphson in several dimensions
%-----------

% calculate first mismatch

R1 = FeT - Feprime - FeL1 - FeL2 - FeL3;
R2 = K1 * L1 * Feprime - FeL1;
R3 = K2 * L2 * Feprime - FeL2;
R4 = K3 * L3 * Feprime - FeL3;
rhs = [R1;R2;R3;R4];
err0 = sum( rhs.^2 );
fprintf('Step 0: Error %e\n',err0)

% do one N-R step and check error

jac = [-1, -1, -1, -1;...
      K1*L1, -1 - K1*Feprime, 0, 0;...
      K2*L2, 0, -1 - K2*Feprime, 0;...
      K3*L3, 0, 0, -1 - K3*Feprime];

dsol = -jac\rhs;
sol = sol0 + dsol;

Feprime = sol(1);
FeL1    = sol(2); 
FeL2    = sol(3);
FeL3    = sol(4); 

L1 = L1T - FeL1;
L2 = L2T - FeL2;
L3 = L3T - FeL3;

R1 = FeT - Feprime - FeL1 - FeL2- FeL3;
R2 = K1 * L1 * Feprime - FeL1;
R3 = K2 * L2 * Feprime - FeL2;
R4 = K3 * L3 * Feprime - FeL3;
rhs = [R1;R2;R3;R4];
err = sum( rhs.^2 );
fprintf('Step 1: Error %e\n',err)

for n=2:10;
  jac = [-1, -1, -1, -1;...
      K1*L1, -1 - K1*Feprime, 0, 0;...
      K2*L2, 0, -1 - K2*Feprime, 0;...
      K3*L3, 0, 0, -1 - K3*Feprime];

  dsol = -jac\rhs;
  sol = sol + dsol;
  Feprime = sol(1);
  FeL1    = sol(2); 
  FeL2    = sol(3);
  FeL3    = sol(4); 
  
  L1 = L1T - FeL1;
  L2 = L2T - FeL2;
  L3 = L3T - FeL3;

  R1 = FeT - Feprime - FeL1 - FeL2- FeL3;
  R2 = K1 * L1 * Feprime - FeL1;
  R3 = K2 * L2 * Feprime - FeL2;
  R4 = K3 * L3 * Feprime - FeL3;
  rhs = [R1;R2;R3;R4];
  err = sum( rhs.^2 );
  fprintf('Step %i: Error %e\n',n,err)

end

% print speciation and mass balances:
fprintf('Fe'':  %8.4f\n',Feprime)
fprintf('FeL1: %8.4f\n',FeL1)
fprintf('FeL2: %8.4f\n',FeL2)
fprintf('FeL3: %8.4f\n',FeL3)

% check whether speciation fits with the original equations
R1 = FeT - Feprime - FeL1 - FeL2;
R2 = K1 * L1 * Feprime - FeL1;
R3 = K2 * L2 * Feprime - FeL2;

