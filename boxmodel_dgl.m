function dydt  = boxmodel_dgl(y,~)
% function dydt  = boxmodel_dgl(y,~)
%
% definition of the rhs of the differential equations of a 14-box
% global box model for anayzing feedbacks in marine iron and ligand
% concentrations

global params

% 12 by 12 matrix of the right hand side from advection
advect = params.advect; 
% ligand stability constant
KFeL = params.KFeL;
% ligand degradation rate
lambda = params.lambda;

nbox = 12;
fe   = y(1,1:12);
lig  = y(1,13:24);
po4  = y(1,25:36);

% first calculate the rate of change from advection
dfedt  = advect*fe;
dligdt = advect*lig;
dpo4dt = advect*po4;

% for ligands, there is a uniform degradation rate
dligdt = dligdt -lambda*lig;

% calculate free iron

% calculate biological productivity in the 7 surface layers

% in the end, glue all rates of change into one long vector
dydt = [dfedt; dligdt; dpo4dt]

return
