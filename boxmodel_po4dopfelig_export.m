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
lig_init = 1.0 + zeros(size(po4_init));
conc_init = [po4_init;dop_init;fe_init;lig_init];

% integrate
tspan = (0:50:3000);
conc = ode23s(@boxmodel_dgl_po4dopfelig_export, tspan, conc_init);

% plots of time development

figure(1)
plot(conc.x,conc.y(1:12,:));
figure(2)
plot(conc.x,conc.y(13:24,:));
figure(3)
plot(conc.x,conc.y(25:36,:));
figure(4)
plot(conc.x,conc.y(37:48,:));

% some diagnostics
totpo4 = params.volume' * conc.y(1:12,:);
totfe = params.volume' * conc.y(25:36,:);
c_export = export * params.redfield_c2p * 12 * 1.0e-18;

for k=1:12
  fprintf('Box %i: %s\n',k,params.long_names{k});
  fprintf('  DFe: %5.2f, Fe: %5.3f\n',conc.y(k+24,end), feprime(k));
  fprintf('  PO4: %5.2f, obs. PO4: %5.2f\n',conc.y(k,end),params.po4init(k));
  fprintf('  DOP: %5.2f\n',conc.y(k+12,end));
  fprintf('  Lig: %5.2f\n',conc.y(k+36,end));
  if k<=5,
    fprintf('  Export: %5.2f PgC \n', c_export(k));
  end
  fprintf('\n');
end

do_plot=1;

if (do_plot),
% plot phosphorus vs data

% plot phosphorus for each box, with data
figure
h = plot(conc.y(1:12,end),'kx');
set(h,'MarkerSize',10, 'LineWidth',2);
hold on
hd = plot(params.po4init,'rx');
set(hd,'MarkerSize',10, 'LineWidth',2);
ylabel('PO_4 [\mumol L^{-1}]');
set(gca,'XTick',(1:12),'XTickLabel',params.names,'XTickLabelRotation',45.0);
set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 3.5]);

% read in iron data and make a plot of data with observations and model
figure
sort_fe_data_into_boxes;
hold on
h = plot(conc.y(25:36,end),'kx');
set(h,'MarkerSize',10, 'LineWidth',2);
ylabel('dFe [nmol L^{-1}]');
set(gca,'XTick',(1:12),'XTickLabel',params.names,'XTickLabelRotation',45.0);
set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 2]);
% print('fe_vs_data_po4dopdfe_nosed.png','-dpng');

end
