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

% read box model initial parameters, to have the bix manes and the constant 
% ligand concentration for run 0L
global params
boxmodel_init_params

totallig = [params.lig*ones(12,1), a(:,4), b(:,4)+b(:,5), c(:,4)+c(:,5)];

lig1 = [a(:,4), b(:,4), c(:,4)];
lig2 = [zeros(12,1), b(:,5), c(:,5)];
bothligs = zeros(12,3,2);
bothligs(:,:,1) = lig1;
bothligs(:,:,2) = lig2;

addpath('~/matlab/tools/plotBarStackGroups/')
hb = plotBarStackGroups(bothligs,params.names);
%set(hb(1,1),'EdgeColor','none')
set(hb(1,1),'FaceColor',[0.0 0.15 0.70],'EdgeColor','none')
set(hb(1,2),'EdgeColor','none')
%set(hb(2,1),'EdgeColor','none')
set(hb(2,1),'FaceColor',[0.48 0.75 0.026],'EdgeColor','none')
set(hb(2,2),'FaceColor',[0.30 0.50 0.016],'EdgeColor','none')
%set(hb(3,1),'EdgeColor','none')
set(hb(3,1),'FaceColor',[0.96 0.80 0.37],'EdgeColor','none')
set(hb(3,2),'FaceColor',[0.64 0.52 0.24],'EdgeColor','none')
hold on
hl = plot([0 13],[params.lig params.lig],'k--');
ylabel('Total ligand [nmol L^{-1}]');
set(gca,'XTickLabelRotation',45.0);
set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 1.8],'box','on');

print('ligands_allruns.png','-dpng','-r600');
