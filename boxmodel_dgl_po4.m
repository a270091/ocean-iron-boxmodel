function dydt  = boxmodel_dgl_po4(~,y)
% function dydt  = boxmodel_dgl(y,~)
%
% definition of the rhs of the differential equations of a 12-box
% global box model for anayzing feedbacks in marine iron and ligand
% concentrations

global params

nbox = 12;
po4  = y(1:12,1);

% 12 by 12 matrix of the right hand side from advection
volume = params.volume;
area   = params.area;
advect = params.advect; 
export = params.export;
remin  = (params.reminfrac * export) ./ volume;
uptake = zeros(12,1);
uptake(1:5) = export ./ volume(1:5);

% calculate the rate of change from advection, biological uptake and
% remineralization 
dpo4dt = advect*po4 - uptake + remin;

% for ligands, there is a uniform degradation rate
% dligdt = dligdt -lambda*lig;

% calculate free iron

% calculate biological productivity in the 7 surface layers

% in the end, glue all rates of change into one long vector
dydt = [dpo4dt];

return
