function dydt  = boxmodel_dgl_po4fe(~,y)
% function dydt  = boxmodel_dgl(~,y)
%
% definition of the rhs of the differential equations of a 12-box
% global box model for anayzing feedbacks in marine iron and ligand
% concentrations
% this is the verion for computing concentrations of PO4 and Fe, 
% assuming a constant uniform ligand concentration

global params

global rhs feprime

nbox = 12;
po4  = y(1:nbox,1);
dfe  = y((nbox+1):(2*nbox),1);

% 12 by 12 matrix of the right hand side from advection
volume = params.volume;
area   = params.area;
advect = params.advect; 
export = params.export;
remin  = (params.reminfrac * export) ./ volume;

% iron parameters
rfe2p  = params.rfe2p;
lig    = params.lig;  % for constant ligand concentration
klig   = params.klig; % 
kscav  = params.kscav;
dust   = params.dust;
sol    = params.dustsol;
hydro  = params.hydrothermal;


% calculate free iron 
p = lig - dfe + 1/klig;
q = dfe/klig;
feprime = -p/2 + sqrt(q + (p/2).^2);

uptake = zeros(12,1);
uptake(1:5) = export ./ volume(1:5);

% calculate the rate of change from advection, biological uptake and
% remineralization 
dpo4dt = advect*po4 - uptake + remin;
ddfedt = advect*dfe - rfe2p*uptake + rfe2p*remin - kscav*feprime + ...
    sol*dust./volume + hydro;

% save individual terms on the rhs of the Fe equation for analysis
rhs.advect = advect*dfe;
rhs.uptake = -rfe2p*uptake;
rhs.remin  = rfe2p*remin;
rhs.scav   = -kscav*feprime;
rhs.dust   = sol*dust./volume;
rhs.hydro   = hydro;

% in the end, glue all rates of change into one long vector
dydt = [dpo4dt; ddfedt];

return
