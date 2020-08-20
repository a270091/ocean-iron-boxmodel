function [Feprime,nit] = func_feprime_three_ligands(dfe,l1t,l2t,l3t,iverb)
%--------------------------------------------------------------------------
% function form of the calculation of the ligand equilibrium with
% three different ligands. Ligand strength is obtained via global variables
%--------------------------------------------------------------------------

global params

if (nargin<5), iverb=0; end

FeT = dfe;
L1T = l1t;
L2T = l2t;
L3T = l3t;
K1  = params.klig;
K2  = params.ksid;
K3  = params.khum;

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
if iverb, fprintf('Step 0: Error %e\n',err0); end

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
if iverb, fprintf('Step 1: Error %e\n',err); end

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
  if iverb, fprintf('Step %i: Error %e\n',n,err); end
  
  if (err<1.0e-12), break; end
end

nit = n;

return

