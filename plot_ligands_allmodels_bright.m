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

% read box model initial parameters, to have the bix manes and the constant 
% ligand concentration for run 0L
global params
boxmodel_init_params

% now decide which of the model runs: 
% CL: all ligands == 1 
% 1L: one ligand, case a
% 1H: one ligand, case d2
% 2L: one ligand, case c
% 2H: one ligand, case d


totallig = [params.lig*ones(12,1), a(:,4), d2(:,4), c(:,4)+c(:,5), ...
            d(:,4)+d(:,5)];

lig1 = [a(:,4), d2(:,4), c(:,4),d(:,4)];
lig2 = [zeros(12,1), zeros(12,1), c(:,5), d(:,5)];
bothligs = zeros(12,4,2);
bothligs(:,:,1) = lig1;
bothligs(:,:,2) = lig2;

%------------------------------------------------------
% define colors for the 'Bright' color scheme 
% (https://personal.sron.nl/~pault/)
bpurple = [170, 51,119]/256;
bred    = [238,102,119]/256;
byellow = [204,187, 68]/256;
bgreen  = [ 34,136, 51]/256;
bcyan   = [102,204,238]/256;
bblue   = [ 86,119,170]/256;

% light version of these colors
fac = 0.5;
byellowl = byellow + fac * (1 - byellow); 
bcyanl   = bcyan   + fac * (1 - bcyan); 

addpath('~/matlab/tools/plotBarStackGroups/')
hb = plotBarStackGroups(bothligs,params.names);
%set(hb(1,1),'EdgeColor','none')
set(hb(1,1),'FaceColor',bblue,'EdgeColor','none')
set(hb(1,2),'EdgeColor','none')
set(hb(2,1),'FaceColor',bgreen,'EdgeColor','none')
set(hb(2,2),'EdgeColor','none')
%set(hb(2,1),'EdgeColor','none')
set(hb(3,1),'FaceColor',bcyan,'EdgeColor','none')
set(hb(3,2),'FaceColor',bcyanl,'EdgeColor','none')
%set(hb(3,1),'EdgeColor','none')
set(hb(4,1),'FaceColor',byellow,'EdgeColor','none')
set(hb(4,2),'FaceColor',byellowl,'EdgeColor','none')
hold on
hl = plot([0 13],[params.lig params.lig],'k--');
ylabel('Total ligand [nmol L^{-1}]');
set(gca,'XTickLabelRotation',45.0);
set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 1.8],'box','on');

print('ligands_allruns_may21.png','-dpng','-r600');
