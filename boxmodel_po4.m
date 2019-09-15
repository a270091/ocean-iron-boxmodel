% to test out, first a model for phosphate only

% initialize model parameters
global params
boxmodel_init_params()

% initial PO4 distribution
po4_init = params.po4init;

% integrate
tspan = [0:50:3000];
po4 = ode23(@boxmodel_dgl_po4, tspan, po4_init);

% plot
plot(po4.x,po4.y);

% some diagnostics
totpo4 = params.volume' * po4.y;

