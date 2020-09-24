% optimizes parameters for the ligand dynamics to obtain a better
% fit with iron data

%----
% first step: read iron data and create median values for the boxes
%----
sort_fe_data_into_boxes

global fe_data pvec_dimensional
fe_data = femedian

%----
% initial guess for parameters (for normalization, we define the optimized
% parameters as simple factors multiplying the original dimensional parameters)
%----

pvec_ini = ones(3,1);
pvec_dimensional(1) = 0.3;      % params.hum
pvec_dimensional(2) = 2.5e-4 * 116; % Lig to P ratio in POC remineralization
                                    % (Lig:C ratio * Redfield C:P)
pvec_dimensional(3) = 1.8 * 2.5e-4 * 116 * 0.67;   % params.rlig2p2
pvec = pvec_ini;

[f_ini,dfe_ini] = costf2_hum_boxmodel_po4dopfe3lig_export(pvec_ini);

%----
% search for better parameters. Unlike the other optimization, 
% here we just do a latin hypercube
%----
nk = 11;
nm = 9;
misfit = zeros(nk,nm);
p1 = zeros(nk,1);
p1max = 2.0;
p2 = zeros(nk,1);
for k=1:nk
  pvec(1) = (k-1)*p1max/(nk-1);
  p1(k) = pvec(1);
  for m=1:nm
    fprintf('%i %i\n',k,m)
    pvec(2) = (m-1)/(nm-1)*p1max;
    pvec(3) = pvec(2);
    p2(m) = pvec(2);
    misfit(k,m) = costf2_hum_boxmodel_po4dopfe3lig_export(pvec);
  end
end

figure;
plot(p1,misfit(:,1));
hold on
plot(p1,misfit(:,2),'r');
plot(p1,misfit(:,3),'g');
plot(p1,misfit(:,4),'r--');
plot(p1,misfit(:,5),'g--');
plot(p1,misfit(:,6),'r');
plot(p1,misfit(:,7),'g');
plot(p1,misfit(:,8),'r--');
plot(p1,misfit(:,9),'g--');

return

