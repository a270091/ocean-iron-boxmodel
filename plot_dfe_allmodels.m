% read equilibrium concentrations from the model runs

fid = fopen('results/equil_po4dopfe_export.dat','r');
a0 = zeros(12,3);
for k=1:12
  a0(k,:) = fscanf(fid,'%f',3);
end
fclose(fid);

fid = fopen('results/equil_po4dopfelig_export.dat','r');
a = zeros(12,4);
for k=1:12
  a(k,:) = fscanf(fid,'%f',4);
end
fclose(fid);

fid = fopen('results/equil_po4dopfe2lig_export_2l1.dat','r');
b = zeros(12,5);
for k=1:12
  b(k,:) = fscanf(fid,'%f',5);
end
fclose(fid);

fid = fopen('results/equil_po4dopfe2lig_export_2l2.dat','r');
c = zeros(12,5);
for k=1:12
  c(k,:) = fscanf(fid,'%f',5);
end
fclose(fid);

% extract PO4 from the model runs
dfe_cl  = a0(:,3);
dfe_1l  = a(:,3);
dfe_2l1 = b(:,3);
dfe_2l2 = c(:,3);

% read box model initial parameters, to have the box names and the 
% PO4 distribution from WOA
global params
boxmodel_init_params

figure(1)
clf
sort_fe_data_into_boxes;

hold on
h0 = plot(dfe_cl,'.');
set(h0,'MarkerSize',30, 'LineWidth',2,'Color',[0.8 0.2 0.1]);
h1 = plot(dfe_1l,'.');
set(h1,'MarkerSize',40, 'LineWidth',2,'Color',[0 0.15 0.70]);
h2 = plot(dfe_2l1,'.');
set(h2,'MarkerSize',30, 'LineWidth',2,'Color',[0.48 0.75 0.26]);
h3 = plot(dfe_2l2,'.');
set(h3,'MarkerSize',20, 'LineWidth',2,'Color',[0.960 0.80 0.37]);

ylabel('dFe [nmol L^{-1}]');
set(gca,'XTick',(1:12),'XTickLabel',params.names,'XTickLabelRotation',45.0);
set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 2]);

print('fe_vs_data_allmodels.png','-dpng','-r600');
return

