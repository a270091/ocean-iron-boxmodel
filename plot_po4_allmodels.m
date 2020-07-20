% read equilibrium concentrations from the model runs

fid = fopen('equil_po4dopfe_export.dat','r');
a0 = zeros(12,3);
for k=1:12
  a0(k,:) = fscanf(fid,'%f',3);
end
fclose(fid);

fid = fopen('equil_po4dopfelig_export.dat','r');
a = zeros(12,4);
for k=1:12
  a(k,:) = fscanf(fid,'%f',4);
end
fclose(fid);

fid = fopen('equil_po4dopfe2lig_export_2l1.dat','r');
b = zeros(12,5);
for k=1:12
  b(k,:) = fscanf(fid,'%f',5);
end
fclose(fid);

fid = fopen('equil_po4dopfe2lig_export_2l2.dat','r');
c = zeros(12,5);
for k=1:12
  c(k,:) = fscanf(fid,'%f',5);
end
fclose(fid);

% extract PO4 from the model runs
po4_cl  = a0(:,1);
po4_1l  = a(:,1);
po4_2l1 = b(:,1);
po4_2l2 = c(:,1);

% read box model initial parameters, to have the box names and the 
% PO4 distribution from WOA
global params
boxmodel_init_params
po4_woa = params.po4init;

figure(1)
h0 = plot(po4_woa, po4_cl,'x');
set(h0,'LineWidth',4,'MarkerSize',10,'Color',[0.6 0.2 0.1]);
hold on 
h1 = plot(po4_woa, po4_1l,'x');
set(h1,'LineWidth',3,'MarkerSize',10,'Color',[0 0.4470 0.7410]);
h2 = plot(po4_woa, po4_2l1,'x');
set(h2,'LineWidth',3,'MarkerSize',10,'Color',[0.8290 0.5940 0.0250]);
h3 = plot(po4_woa, po4_2l2,'x');
set(h3,'LineWidth',3,'MarkerSize',10,'Color',[0.3660 0.5740 0.0880]);
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',12);
xlabel('WOA-derived PO_4');
ylabel('model PO_4');
yl = get(gca,'YLim');
hd = plot(yl,yl,'k');

