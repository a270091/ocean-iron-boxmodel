function [f, po4_init, po4_final] = costf_exp_boxmodel_po4(pvec,rplot)
% f = costf_boxmodel_po4(pvec)

% initialize model parameters
global params
boxmodel_init_params()

if (nargin>0), 
  redfield_c2p = 106.0;
  params.export = pvec(1:5) * 1.0e18/12/redfield_c2p;
end

% initial PO4 distribution
po4_init = params.po4init;

% integrate
tspan = [0:50:3000];
po4 = ode23(@boxmodel_dgl_po4, tspan, po4_init);

% calculate RMS difference between final and initial PO4
po4_final = po4.y(:,end);
diff = po4_init - po4_final;
f = sqrt( sum( diff.^2 ) );

if (nargin>1),
    plot(po4.x,po4.y);
    po4.y(:,end);
end

return
end
