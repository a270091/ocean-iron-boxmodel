function dydt  = boxmodel_dgl_po4dopfelig_export(~,y)
% function dydt  = boxmodel_dgl(~,y)
%
% definition of the rhs of the differential equations of a 12-box
% global box model for anayzing feedbacks in marine iron and ligand
% concentrations
% this is the version for computing concentrations of PO4, DOP, Fe, 
% and a spatially varying ligand concentration

global params

global rhs feprime export

nbox = 12;
po4  = y(1:nbox,1);
dop  = y((nbox+1):(2*nbox),1);
dfe  = y((2*nbox+1):(3*nbox),1);
lig  = y((3*nbox+1):(4*nbox),1);

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

% ligand parameters
rlig2p  = params.rlig2p;
rlig2p2 = params.rlig2p2;
ligrem  = params.ligrem;
surffac = params.ligfac;

% calculate free iron 
p = lig - dfe + 1/klig;
q = dfe/klig;
feprime = -p/2 + sqrt(q + (p/2).^2);

% calculate the rate of change from advection, biological uptake and
% remineralization 
dpo4dt = advect*po4 - uptake + remin + dopremin*dop;
ddopdt = advect*dop + dopfrac*uptake - dopremin*dop;
ddfedt = advect*dfe - rfe2p*uptake + rfe2p*remin  ...
    + rfe2p*dopremin*dop - kscav*feprime + sol*dust./volume + ...
	  hydfac*hydro + sedfac*sedfe;
ligremin = ligrem*lig;
ligremin(1:5) = ligremin(1:5)*surffac;
dligdt = advect*lig + rlig2p*remin + rlig2p2*uptake - ligremin; 

% save individual terms on the rhs of the Fe equation for analysis
rhs.advect = advect*dfe;
rhs.uptake = -rfe2p*uptake;
rhs.remin  = rfe2p*remin;
rhs.scav   = -kscav*feprime;
rhs.dust   = sol*dust./volume;
rhs.hydro  = hydfac*hydro;
rhs.sedfe  = sedfac*sedfe;
rhs.lig1   = rlig2p*remin;
rhs.lig2   = rlig2p2*uptake;
rhs.lig3   = -ligremin;

% in the end, glue all rates of change into one long vector
dydt = [dpo4dt; ddopdt; ddfedt; dligdt];

return

end

