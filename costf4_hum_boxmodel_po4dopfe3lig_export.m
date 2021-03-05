function [f,dfe_final] = costf4_hum_boxmodel_po4dopfe3lig_export(pvec)

global fe_data pvec_dimensional

% global fields for analysis
global params
global rhs feprime export

% initialize model parameters
boxmodel_init_params()

% adapt some parameters for ligands
if (nargin>0) 
  params.hum      = pvec_dimensional(1)*pvec(1);
  params.rlig2p   = pvec_dimensional(2)*pvec(2);
  params.ligfac   = pvec_dimensional(3)*pvec(3);
  params.ligrem   = pvec_dimensional(4)*pvec(4);
end
params.sidremin = 0.5*0.1;
params.KFe_bact = 0.1*0.01;

% initial PO4 distribution
po4_init = params.po4init;
dop_init = zeros(size(po4_init));
fe_init = 0.6 + zeros(size(po4_init));
lig_init = 1.0 + zeros(size(po4_init));
sid_init = zeros(size(po4_init));
conc_init = [po4_init;dop_init;fe_init;lig_init;sid_init];

% integrate
tspan = (0:50:5000);
conc = ode23s(@boxmodel_dgl_po4dopfe3lig_export, tspan, conc_init);

% calculate RMS difference between measured and modeled Fe
dfe_final = conc.y(25:36,end);
diff = fe_data - dfe_final;
% f = sqrt( sum( (diff.^2).*params.volume ) ./ sum(params.volume) );
f = sqrt( sum( (diff.^2) ) );
%f = sum( abs(diff ) );

if min(pvec)<0
    f = 10+sum(abs(pvec));
end

return
end
