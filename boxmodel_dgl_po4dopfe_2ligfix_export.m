function dydt  = boxmodel_dgl_po4dopfe_export(~,y)
% function dydt  = boxmodel_dgl(~,y)
%
% definition of the rhs of the differential equations of a 12-box
% global box model for anayzing feedbacks in marine iron and ligand
% concentrations
% this is the verion for computing concentrations of PO4 and Fe, 
% assuming a constant uniform ligand concentration

global params
global rhs feprime export
global ligfix sidfix

nbox = 12;
po4  = y(1:nbox,1);
dop  = y((nbox+1):(2*nbox),1);
dfe  = y((2*nbox+1):(3*nbox),1);
lig = ligfix;
sid = sidfix;

% 12 by 12 matrix of the right hand side from advection
volume = params.volume;
area   = params.area;
advect = params.advect; 

% dop parameters
dopfrac  = params.dopfrac;
dopremin = params.dopremin;

% prognostic calculation of export
uptake = zeros(12,1);
uptake(1:5) = params.mumax .* po4(1:5) .* (dfe(1:5) ./ ...
    (dfe(1:5) + params.kFe)); 
export = uptake(1:5)* (1 - dopfrac) .* volume(1:5);
remin  = (params.reminfrac * export) ./ volume;

% iron parameters
rfe2p  = params.rfe2p;
klig   = params.klig; % 
kscav  = params.kscav;
dust   = params.dust;
sol    = params.dustsol;
hydro  = params.hydrothermal;
hydfac = params.hydro_fac;
sedfe  = params.sediment_fe;
sedfac = params.sed_fac;
ksid = params.ksid; 

% calculate free iron with two different ligands from third-order
% polynomial equation
feprime = zeros(12,1);
for k=1:nbox
  % polynomial coefficients
  a3 = -klig*ksid;
  a2 = ( klig*ksid*(dfe(k) - lig(k) - sid(k)) - klig - ksid );
  a1 = (klig + ksid)*dfe(k) - 1 - klig*lig(k) - ksid*sid(k);
  a0 = dfe(k);
  % solve cubic equation and return the largest real root
  test = sort( roots([a3 a2 a1 a0]), 'ComparisonMethod','real' );
  % calculate speciation
  feprime(k) = test(end);
end
% fprintf('feprime: %i %i\n', size(feprime))

%uptake = zeros(12,1);
%uptake(1:5) = export ./ volume(1:5);

% calculate the rate of change from advection, biological uptake and
% remineralization 
dpo4dt = advect*po4 - uptake + remin + dopremin*dop;
ddopdt = advect*dop + dopfrac*uptake - dopremin*dop;
ddfedt = advect*dfe - rfe2p*uptake + rfe2p*remin  ...
    + rfe2p*dopremin*dop - kscav*feprime + sol*dust./volume + ...
	  hydfac*hydro + sedfac*sedfe;

% save individual terms on the rhs of the Fe equation for analysis
rhs.advect = advect*dfe;
rhs.uptake = -rfe2p*uptake;
rhs.remin  = rfe2p*remin;
rhs.scav   = -kscav*feprime;
rhs.dust   = sol*dust./volume;
rhs.hydro  = hydfac*hydro;
rhs.sedfe  = sedfac*sedfe;

% in the end, glue all rates of change into one long vector
dydt = [dpo4dt; ddopdt; ddfedt];

return
