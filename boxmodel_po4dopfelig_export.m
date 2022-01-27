% to test out, first a model for phosphate only

tic;
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
tspan = (0:50:5000);
conc = ode23s(@boxmodel_dgl_po4dopfelig_export, tspan, conc_init);

toc;

% plots of time development

% figure(1)
% plot(conc.x,conc.y(1:12,:));
% figure(2)
% plot(conc.x,conc.y(13:24,:));
% figure(3)
% plot(conc.x,conc.y(25:36,:));
% figure(4)
% plot(conc.x,conc.y(37:48,:));

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
    fprintf('  Export: %5.2f PgC/yr \n', c_export(k));
  end
  fprintf('\n');
end

% calculate some global numbers and print them

dfe_av = sum( conc.y(25:36,end) .* params.volume ) / sum(params.volume);
dfe_tot = sum( conc.y(25:36,end) .* params.volume ) * 1.0e-6; 
dust   = params.dust;
sol    = params.dustsol;
dust_fe_sol = sum( sol*dust ) * 1.0e-6;
sedfe  = params.sediment_fe;
sedfac = params.sed_fac;
sed_fe_sol = sum(sedfac * sedfe .* params.volume ) * 1.0e-6;
hydro  = params.hydrothermal;
hydfac = params.hydro_fac;
hyd_fe_sol = sum(hydfac * hydro .* params.volume ) * 1.0e-6;
residence_time = dfe_tot / (dust_fe_sol + sed_fe_sol + hyd_fe_sol);
fprintf('Total export: %5.2f PgC/yr \n', sum(c_export));
fprintf('Average dFe concentration: %5.2f nmol/L \n',dfe_av);
fprintf('Dust iron input: %8.2e mol/yr \n',dust_fe_sol);
fprintf('Sediment iron input: %8.2e mol/yr \n',sed_fe_sol);
fprintf('Hydrothermal iron input: %8.2e mol/yr \n',hyd_fe_sol);
fprintf('Residence time: %6.2f yr \n',residence_time);

% save final concentrations as a table
% (for later plotting of equilibrium ligands)
finalstate = conc.y(:,end);
fid = fopen('equil_po4dopfelig_export.dat','w');
for k=1:12
  fprintf(fid,'%8.4f %8.3f %8.4f %8.4f\n', finalstate(k),...
	 finalstate(k+12),finalstate(k+24),finalstate(k+36));
end
fclose(fid);

% save parameter values as a matlab-file
save('parameters_1l.mat','-struct','params');

do_plot=0;

if (do_plot),
% plot phosphorus vs data
figure(1)
po4_end = conc.y(1:12,end);
hh = plot(po4_init, po4_end,'x');
hold on 
set(hh,'LineWidth',2,'MarkerSize',12,'Color',[0 0.6 0.3]);
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',12);
xlabel('WOA-derived PO_4');
ylabel('model PO_4');
yl = get(gca,'YLim');
hd = plot(yl,yl,'k');
print('PO4_boxmodel_po4dopfelig_sed.png','-dpng')

% plot phosphorus for each box, with data
%figure
%h = plot(conc.y(1:12,end),'kx');
%set(h,'MarkerSize',10, 'LineWidth',2);
%hold on
% hd = plot(params.po4init,'rx');
% set(hd,'MarkerSize',10, 'LineWidth',2);
% ylabel('PO_4 [\mumol L^{-1}]');
% set(gca,'XTick',(1:12),'XTickLabel',params.names,'XTickLabelRotation',45.0);
% set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 3.5]);

% read in iron data and make a plot of data with observations and model
figure(2)
plot_geotraces_idp2;
hold on
h = plot(conc.y(25:36,end),'kx');
set(h,'MarkerSize',10, 'LineWidth',2);
ylabel('dFe [nmol L^{-1}]');
set(gca,'XTick',(1:12),'XTickLabel',params.names,'XTickLabelRotation',45.0);
set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 2]);
print('fe_vs_data_po4dopdfelig_sed.png','-dpng');

end
