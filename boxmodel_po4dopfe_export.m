% to test out, first a model for phosphate only

% global fields for analysis
global params
global rhs feprime export

% initialize model parameters
boxmodel_init_params()

% initial PO4 distribution
po4_init = params.po4init;
dop_init = zeros(size(po4_init));
fe_init = 0.6 + zeros(size(po4_init));
conc_init = [po4_init;dop_init;fe_init];

% integrate
tspan = [0:50:3000];
conc = ode23(@boxmodel_dgl_po4dopfe_export, tspan, conc_init);

% plot
plot(conc.x,conc.y(1:12,:));
figure
plot(conc.x,conc.y(13:24,:));
figure
plot(conc.x,conc.y(25:36,:));

% some diagnostics
totpo4 = params.volume' * conc.y(1:12,:);
totfe = params.volume' * conc.y(25:36,:);

for k=1:12
  fprintf('Box %i: %s\n',k,params.long_names{k});
  fprintf('  DFe: %5.2f, Fe: %5.3f\n',conc.y(k+24,end), feprime(k));
  fprintf('\n');
end


