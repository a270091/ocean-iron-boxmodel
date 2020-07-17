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

% read box model initial parameters, to have the bix manes and the constant 
% ligand concentration for run 0L
global params
boxmodel_init_params

totallig = [params.lig*ones(12,1), a(:,4), b(:,4)+b(:,5), c(:,4)+c(:,5)];

bar(totallig)
%set(h,'MarkerSize',10, 'LineWidth',2);
ylabel('Total ligand [nmol L^{-1}]');
set(gca,'XTick',(1:12),'XTickLabel',params.names,'XTickLabelRotation',45.0);
set(gca,'FontSize',12,'XLim',[0.5,12.5],'YLim',[0 2]);
% print('fe_vs_data_po4dopdfe2lig_sed.png','-dpng');
