function boxmodel_init_params()

global params

% names for the boxes
params.long_names = {'North Pacific Surface','Pacific Equatorial Surface',...
		    'Southern Ocean Surface','Atlantic Equatorial Surface',...
		    'North Atlantic Surface','Pacific Intermediate Water',...
		    'Atlantic Intermediate Water','Pacific Deep Water',...
		    'Deep Southern Ocean','North Atlantic Deep Water',...
		    'Pacific Bottom Water','Atlantic Bottom Water'};
params.names = {'NPS','EQPS','SOS','EQAS','NAS','PIW','AIW',...
		'PDW','SOD','NADW','PBW','ABW'};

% volumes for the boxes in m^3, and surface area for the 
% upper 5 boxes in m^2

params.volume = zeros(12,1);
params.volume(1) = 3.817e+15;
params.volume(2) = 1.562e+16;
params.volume(3) = 1.155e+16;
params.volume(4) = 4.276e+15;
params.volume(5) = 2.578e+15;
params.volume(6) = 1.789e+17;
params.volume(7) = 4.486e+16;
params.volume(8) = 5.486e+17;
params.volume(9) = 7.830e+16;
params.volume(10) = 1.925e+17;
params.volume(11) = 2.661e+17;
params.volume(12) = 9.946e+16;

params.area = zeros(5,1);
params.area(1) = 3.444e+13;
params.area(2) = 1.404e+14;
params.area(3) = 1.036e+14;
params.area(4) = 3.850e+13;
params.area(5) = 2.324e+13;

params.po4init = zeros(12,1);
params.po4init(1)  = 0.73;
params.po4init(2)  = 0.43;
params.po4init(3)  = 1.16;
params.po4init(4)  = 0.34;
params.po4init(5)  = 0.32;
params.po4init(6)  = 1.93;
params.po4init(7)  = 1.43;
params.po4init(8)  = 2.64;
params.po4init(9)  = 2.28;
params.po4init(10) = 1.58;
params.po4init(11) = 2.36;
params.po4init(12) = 1.88;

%------------------------------------------------------------------------
% advection: first define a three-column matrix where the first 
% column gives the number of the box that the advection originates
% from, the second the box it goes into, and the third the flux
% in Sverdrup (10^6 m^3/s)
%------------------------------------------------------------------------

adv_matrix = zeros(21,3);
x = 3;
adv_matrix(1,:) =  [9 , 12, 4 ]; % AABW formation in tlantic
adv_matrix(2,:) =  [10, 12, 3 ]; % Entrainment of NADW into AABW
adv_matrix(3,:) =  [12, 10, 7 ]; % AABW mixing into NADW
adv_matrix(4,:) =  [10,  9, 13]; % NADW upwelling directly into deep water formation
adv_matrix(5,:) =  [10,  7, 5 ]; % NADW mixing into AAIW + CW (to balance
                                 % northward surface flow of 18Sv) 
adv_matrix(6,:) =  [ 5, 10, 19]; % NADW formation
adv_matrix(7,:) =  [ 4,  5, 18]; % northward surface flow at 24N
adv_matrix(8,:) =  [ 3,  7,  2];
adv_matrix(9,:) =  [ 7,  4, 18];
adv_matrix(10,:) = [10,  8,  5];
adv_matrix(11,:) = [ 6,  7, 11];
adv_matrix(12,:) = [ 1,  5,  1];
adv_matrix(13,:) = [ 9, 11, 25];
adv_matrix(14,:) = [11,  8, 25];
adv_matrix(15,:) = [ 8,  6,  6];
adv_matrix(16,:) = [ 8,  3, 24];
adv_matrix(17,:) = [ 3,  9, 16];
adv_matrix(18,:) = [ 3,  6,  6];
adv_matrix(19,:) = [ 6,  2,1+x];
adv_matrix(20,:) = [ 2,  1,1+x];
adv_matrix(21,:) = [ 1,  6,  x];

%------------------------------------------------------------------------
% we also need some mixing between the surface boxes and the boxes below; 
% mixing is conceptualized as biderictional flux from box A to box B and
% vice versa.
%------------------------------------------------------------------------
mix_matrix = zeros(7,3);
mix_matrix(1,:) = [6, 1, 31.59];  % 40
mix_matrix(2,:) = [6, 2, 38.13];  % 40
mix_matrix(3,:) = [7, 4, 27.04];  % 17
mix_matrix(4,:) = [10, 5, 12.03]; % 10
mix_matrix(5,:) = [8, 3, 4.12];   % 5
mix_matrix(6,:) = [9, 3, 0.0];    % 20
mix_matrix(7,:) = [10, 3, 23.85]; % 40

tr = zeros(12,12);
for k=1:21,
  i = adv_matrix(k,1);
  j = adv_matrix(k,2);
  tr(i,i) = tr(i,i) - adv_matrix(k,3);
  tr(j,i) = tr(j,i) + adv_matrix(k,3);
end

mix = zeros(12,12);
for k=1:7,
  i = mix_matrix(k,1);
  j = mix_matrix(k,2);
  mix(i,i) = mix(i,i) - mix_matrix(k,3);
  mix(j,j) = mix(j,j) - mix_matrix(k,3);
  mix(j,i) = mix(j,i) + mix_matrix(k,3);
  mix(i,j) = mix(i,j) + mix_matrix(k,3);
end
% tr = tr - tr';

params.advect = zeros(12,12);
unitfac = 1.0e6 * 360 * 86400; % to convert volume transport from 10^6 m^3/s
                               % to m^3/yr
for k=1:12,
  params.advect(k,:) = unitfac * (tr(k,:) + mix(k,:))/ params.volume(k);
end

%------------------------------------------------------------------------
% the net biological changes of PO4 from export / remineralization are
% calculated here; this part of the code will later move into the
% differential equation, as soon as export becomes dependent on PO4 and Fe
%------------------------------------------------------------------------

% PO$ concentration changes, calculated from Export Production in each of the 
% surface boxes

params.redfield_c2p = 106;
params.c_export = zeros(5,1);
% export in PgC/yr, derived from Laws et al. 
params.c_export(1) = 0.99; % 1.87;
params.c_export(2) = 1.88; % 2.83;
params.c_export(3) = 2.38; % 2.66;
params.c_export(4) = 0.72; % 1.20;
params.c_export(5) = 1.44; % 1.53;
% conversion to mmol P/yr
params.export = params.c_export*1.0e18/12/params.redfield_c2p;

%------------------------------------------------------------------------
% parameters affecting DOP
%------------------------------------------------------------------------
params.dopfrac = 0.67;   % Parekh: 0.67
params.dopremin = 0.5;   % Parekh: 0.5/yr

%------------------------------------------------------------------------
% here are now the parameters for making export production
% prognostically dependent on PO4 and iron. Standard parameters are 
% from Parekh et al., 2004:
%------------------------------------------------------------------------
% 
% export(PO4) = mumax * PO4 * Fe / (Ke + kFe)
params.mumax = 12; % 1/month = 12/yr
params.kFe   = 0.2;

%------------------------------------------------------------------------
% reminfrac(i,j) says which fraction of export from box j is remineralized
% in box i (all fractions are between zero and, and the sum of all
% fractions in a column is one)
%------------------------------------------------------------------------
params.reminfrac = zeros(12,5);
params.reminfrac(6,1) = 0.85;
params.reminfrac(8,1) = 0.13;
params.reminfrac(11,1)= 0.02;

params.reminfrac(6,2) = 0.85;
params.reminfrac(8,2) = 0.13;
params.reminfrac(11,2)= 0.02;

params.reminfrac(8,3) = 0.54;
params.reminfrac(9,3) = 0.06;
params.reminfrac(10,3)= 0.36;
params.reminfrac(11,3)= 0.02;
params.reminfrac(12,3)= 0.02;

params.reminfrac(7,4) = 0.85;
params.reminfrac(10,4)= 0.13;
params.reminfrac(12,4)= 0.02;

params.reminfrac(10,5)= 0.97;
params.reminfrac(12,5)= 0.03;

%------------------------------------------------------------------------
% now parameters that are relevant for the iron cycle
%------------------------------------------------------------------------
params.rfe2p  = 0.53;    % unit: nmol Fe/umol P, or mmol Fe/mol P, 
                         % calculated from Fe:C = 5.0e-3 mmol/mol and C:P = 106 mol/mol 
params.lig    = 0.6;     % constant ligand concentration in nmol/L
params.klig   = 200.0;   % unit: L/nmol. Corresponds to a stability constant of 2.0e11 L/mol
params.kscav  = 1/50;   % unit: 1/year

% dust values used here were derived from Mahowald/Luo dust deposition
% fields, assuming 3.5 weight percent of iron in dust. This is the total
% iron flux, solubility has to be taken into account separately. 
params.dust   = zeros(12,1);
params.dust(1) = 7.083715e+15;
params.dust(2) = 5.312438e+16;
params.dust(3) = 1.914294e+16;
params.dust(4) = 6.325368e+16;
params.dust(5) = 1.063543e+16;

params.dustsol = 0.01;

% hydrothermal Fe flux, already converted into a volumentric concentration
% rate of change (in nmol/L/yr). Values come from conversion of
% hydrothermal flux of 3He, so that the total flux of Helium-3 into the ocean
% is 1000 mol/yr. Conversion is done so that the total hydrothermal Fe flux
% is 9.0e8 mol/yr. Since some minor basins are missing in the box model,
% the total flux is 8.47e8 mol/yr.

params.hydrothermal = zeros(12,1);
params.hydrothermal(1) = 0.00e+00;
params.hydrothermal(2) = 2.89e-06;
params.hydrothermal(3) = 0.00e+00;
params.hydrothermal(4) = 0.00e+00;
params.hydrothermal(5) = 2.14e-04;
params.hydrothermal(6) = 2.08e-04;
params.hydrothermal(7) = 0.00e+00;
params.hydrothermal(8) = 7.62e-04;
params.hydrothermal(9) = 7.87e-04;
params.hydrothermal(10) = 4.21e-04;
params.hydrothermal(11) = 7.75e-04;
params.hydrothermal(12) = 4.36e-04;

% sedimentary  iron flux, already converted into a volumentric concentration
% rate of change (in nmol/L/yr). Values have been estimated from the sediment
% area contained in each of the WOA 1-degree grid boxes, multiplied by a depth-
% dependent Fe flux per square meter and year, with a maximum flux of 
% 2 micromol/m^2/day. This parameterization stems from Aumont et al., 2015.

params.sediment_fe = zeros(12,1);
params.sediment_fe(1) = 2.96e-01;
params.sediment_fe(2) = 1.79e-01;
params.sediment_fe(3) = 6.25e-02;
params.sediment_fe(4) = 1.53e-01;
params.sediment_fe(5) = 3.74e-01;
params.sediment_fe(6) = 1.62e-02;
params.sediment_fe(7) = 2.75e-02;
params.sediment_fe(8) = 2.70e-03;
params.sediment_fe(9) = 1.71e-02;
params.sediment_fe(10) = 5.84e-03;
params.sediment_fe(11) = 1.97e-03;
params.sediment_fe(12) = 2.33e-03;

return


