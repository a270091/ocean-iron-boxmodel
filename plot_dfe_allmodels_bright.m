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

fid = fopen('results/equil_po4dopfe3lig_export_3l1.dat','r');
d = zeros(12,5);
for k=1:12
  d(k,:) = fscanf(fid,'%f',5);
end
fclose(fid);

fid = fopen('results/equil_po4dopfe3lig_export_3l2.dat','r');
d2 = zeros(12,5);
for k=1:12
  d2(k,:) = fscanf(fid,'%f',5);
end
fclose(fid);

% extract iron (3rd column) from the model runs
dfe_cl  = a0(:,3);
dfe_1l  = a(:,3);
dfe_2l1 = b(:,3);
dfe_2l2 = c(:,3);
dfe_3l1 = d(:,3);
dfe_3l2 = d2(:,3);

% read box model initial parameters, to have the box names and the 
% PO4 distribution from WOA
global params
boxmodel_init_params

figure(1)
clf
sort_fe_data_into_boxes;



% whitebg([0.9 0.9 0.9]);

%------------------------------------------------------
% define colors for the 'Bright' color scheme 
% (https://personal.sron.nl/~pault/)
bpurple = [170, 51,119]/256;
bred    = [238,102,119]/256;
byellow = [204,187, 68]/256;
bgreen  = [ 34,136, 51]/256;
bcyan   = [102,204,238]/256;
bblue   = [ 86,119,170]/256;

%------------------------------------------------------
% make plot 
hold on
h0 = plot(dfe_cl,'.');
set(h0,'MarkerSize',40, 'LineWidth',2,'Color',bred);
h1 = plot(dfe_1l,'.');
set(h1,'MarkerSize',50, 'LineWidth',2,'Color',bblue);
h5 = plot(dfe_3l2,'.');
set(h5,'MarkerSize',40, 'LineWidth',2,'Color',bgreen);
%set(h5,'MarkerSize',40, 'LineWidth',2,'Color',[0.70 0.0 0.95]);
h4 = plot(dfe_3l1,'.');
set(h4,'MarkerSize',30, 'LineWidth',2,'Color',byellow);
% set(h4,'MarkerSize',30, 'LineWidth',2,'Color',[0.00 0.75 0.00]);
%h2 = plot(dfe_2l1,'.');
%set(h2,'MarkerSize',30, 'LineWidth',2,'Color',[0.0 0.9 0.9]);
h3 = plot(dfe_2l2,'.');
set(h3,'MarkerSize',20, 'LineWidth',2,'Color',bcyan);
% set(h3,'MarkerSize',20, 'LineWidth',2,'Color',[0.90 0.90 0.0]);

% add a legend
y0 = 1.8;
dy = 0.1;
dx = 0.4;
hl0 = plot(1.5,y0,'.');
ht0 = text(1.5+dx,y0,'CL');
set(hl0,'MarkerSize',40, 'LineWidth',2,'Color',bred);
set(ht0,'FontSize',12)
hl1 = plot(1.5,y0-dy,'.');
ht1 = text(1.5+dx,y0-dy,'1L');
set(hl1,'MarkerSize',50, 'LineWidth',2,'Color',bblue);
set(ht1,'FontSize',12)
hl5 = plot(1.5,y0-2*dy,'.');
ht5 = text(1.5+dx,y0-2*dy,'1H');
%set(hl5,'MarkerSize',40, 'LineWidth',2,'Color',[0.70 0.0 0.95]);
set(hl5,'MarkerSize',40, 'LineWidth',2,'Color',bgreen);
set(ht5,'FontSize',12)
hl4 = plot(1.5,y0-3*dy,'.');
ht4 = text(1.5+dx,y0-3*dy,'2H');
set(hl4,'MarkerSize',30, 'LineWidth',2,'Color',byellow);
set(ht4,'FontSize',12)
% hl2 = plot(1.5,y0-4*dy,'.');
% ht2 = text(1.5+dx,y0-4*dy,'2L1');
% set(hl2,'MarkerSize',30, 'LineWidth',2,'Color',[0.48 0.75 0.26]);
% set(hl2,'MarkerSize',30, 'LineWidth',2,'Color',[0.0 0.9 0.9]);
% set(ht2,'FontSize',12)
hl3 = plot(1.5,y0-4*dy,'.');
ht3 = text(1.5+dx,y0-4*dy,'2L');
% set(hl3,'MarkerSize',20, 'LineWidth',2,'Color',[0.90 0.90 0.00]);
set(hl3,'MarkerSize',20, 'LineWidth',2,'Color',bcyan);
set(ht3,'FontSize',12)

ylabel('dFe [nmol L^{-1}]');
set(gca,'XTick',(1:12),'XTickLabel',params.names,'XTickLabelRotation',45.0);
set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 2]);

print('fe_vs_data_allmodels_may21.png','-dpng','-r600');
return

