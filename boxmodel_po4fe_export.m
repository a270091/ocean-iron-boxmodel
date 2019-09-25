% a model for phosphate and iron with constant ligands and with export
% production prescribed

% global fields for analysis
global params
global rhs feprime export

% initialize model parameters
boxmodel_init_params()

% initial PO4 distribution
po4_init = params.po4init;
fe_init = 0.6 + zeros(size(po4_init));
conc_init = [po4_init;fe_init];

% integrate
tspan = [0:50:3000];
conc = ode23s(@boxmodel_dgl_po4fe_export, tspan, conc_init);

% plot
plot(conc.x,conc.y(1:12,:));
figure
plot(conc.x,conc.y(13:24,:));

% some diagnostics
totpo4 = params.volume' * conc.y(1:12,:);
totfe = params.volume' * conc.y(13:24,:);
c_export = export * params.redfield_c2p * 12 * 1.0e-18;

for k=1:12
  fprintf('Box %i: %s\n',k,params.long_names{k});
  fprintf('  DFe: %5.2f, Fe: %5.3f\n',conc.y(k+12,end), feprime(k));
  fprintf('  PO4: %5.2f\n',conc.y(k,end));
  if k<=5,
    fprintf('  Export: %5.2f PgC \n', c_export(k));
  end
  fprintf('\n');
end

