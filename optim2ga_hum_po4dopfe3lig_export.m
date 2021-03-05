%-------------
% optimizes parameters for the ligand dynamics to obtain a better
% fit with iron data
%
% this version uses a genetic algorithm as coded in the GODLIKE 
% toolbox
%-------------
addpath('/Users/cvoelker/matlab/tools/GODLIKE');

%---
% start timer, to get an idea how long optimzations take
%---
tic;

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

pvec_ini = ones(5,1);
pvec_dimensional(1) = 0.3;      % params.hum
pvec_dimensional(2) = 2.5e-4 * 116; % Lig to P ratio in POC remineralization
                                    % (Lig:C ratio * Redfield C:P)
pvec_dimensional(3) = 1.8 * 2.5e-4 * 116 * 0.67;   % params.rlig2p2
pvec_dimensional(4) = 0.5e-3; % params.ligrem
pvec_dimensional(5) = 100.0;  % params.ligfac
pvec = pvec_ini;

[f_ini,dfe_ini] = costf2_hum_boxmodel_po4dopfe3lig_export(pvec_ini);

%-----------------
% constrained combinatorial optimzation, using GODLIKE
%-----------------
parmin = [0.05 0.05 0.05 0.05];
parmax = [2.0  2.0  2.0  2.0];
optvec = set_options('TolFun',1.0e-4,'MaxIters',100,'Display','on');

[paramopt,fopt] = GODLIKE(@costf2_hum_boxmodel_po4dopfe3lig_export,10,parmin,parmax,'GA',optvec);

fopt2vec(nopt) = fopt2;
paramopt2vec(nopt,:) = paramopt2;

%--
% save results into file
%--
isave = 0;
if (isave),
    fid = fopen('opt_hum_new.out','w')
    for k=1:nn
        fprintf(fid,'%7.4f %7.4f %7.4f %7.4f %7.4f\n',cost_all(k),...
                pvec_all(k,1),pvec_all(k,2),...
                pvec_all(k,3),pvec_all(k,4));
    end
    fclose(fid);
end

%--
% end timer and print out
%--
totaltime = toc/3600;
fprintf('runtime %5.2f hours \n',totaltime);

% [f_end,dfe_end] = costf2_hum_boxmodel_po4dopfe3lig_export(pvec);
