function [f,po4_init,po4_final] = costf_mix_boxmodel_po4(pvec,rplot)
% f = costf_boxmodel_po4(pvec)

% initialize model parameters
global params
boxmodel_init_params()

if (nargin>0), 
% advection: first define a three-column matrix where the first 
% column gives the number of the box that the advection originates
% from, the second the box it goes into, and the third the flux
% in Sverdrup (10^6 m^3/s)
  adv_matrix = zeros(21,3);
  x = 3;
  adv_matrix(1,:) =  [9 , 12, 4 ]; % AABW formation in tlantic
  adv_matrix(2,:) =  [10, 12, 3 ]; % Entrainment of NADW into AABW
  adv_matrix(3,:) =  [12, 10, 7 ]; % AABW mixing into NADW
  adv_matrix(4,:) =  [10,  9, 13]; % NADW upwelling directly into deep water formation
  adv_matrix(5,:) =  [10,  7, 5 ]; % NADW mixing into AAIW + CW (to balance
                                 % northward surface flow of 18Sv) 
  adv_matrix(6,:) =  [ 5, 10, 19]; % NADW formation
  adv_matrix(7,:) =  [ 4,  5, 18]; % northward surface flow at 24N
  adv_matrix(8,:) =  [ 3,  7,  2];
  adv_matrix(9,:) =  [ 7,  4, 18];
  adv_matrix(10,:) = [10,  8,  5];
  adv_matrix(11,:) = [ 6,  7, 11];
  adv_matrix(12,:) = [ 1,  5,  1];
  adv_matrix(13,:) = [ 9, 11, 25];
  adv_matrix(14,:) = [11,  8, 25];
  adv_matrix(15,:) = [ 8,  6,  6];
  adv_matrix(16,:) = [ 8,  3, 24];
  adv_matrix(17,:) = [ 3,  9, 16];
  adv_matrix(18,:) = [ 3,  6,  6];
  adv_matrix(19,:) = [ 6,  2,1+x];
  adv_matrix(20,:) = [ 2,  1,1+x];
  adv_matrix(21,:) = [ 1,  6,  x];

% we also need some mixing between the surface boxes and the boxes below; 
% mixing is conceptualized as biderictional flux from box A to box B and
% vice versa. 
  mix_matrix(1,:) = [6, 1, pvec(1)];
  mix_matrix(2,:) = [6, 2, pvec(2)];
  mix_matrix(3,:) = [7, 4, pvec(3)];
  mix_matrix(4,:) = [10, 5, pvec(4)];
  mix_matrix(5,:) = [8, 3, pvec(5)];
  mix_matrix(6,:) = [9, 3, pvec(6)];
  mix_matrix(7,:) = [10, 3, pvec(7)];

  tr = zeros(12,12);
  for k=1:21,
    i = adv_matrix(k,1);
    j = adv_matrix(k,2);
    tr(i,i) = tr(i,i) - adv_matrix(k,3);
    tr(j,i) = tr(j,i) + adv_matrix(k,3);
  end

  mix = zeros(12,12);
  for k=1:7,
    i = mix_matrix(k,1);
    j = mix_matrix(k,2);
    mix(i,i) = mix(i,i) - mix_matrix(k,3);
    mix(j,j) = mix(j,j) - mix_matrix(k,3);
    mix(j,i) = mix(j,i) + mix_matrix(k,3);
    mix(i,j) = mix(i,j) + mix_matrix(k,3);
  end
  % tr = tr - tr';

  params.advect = zeros(12,12);
  unitfac = 1.0e6 * 360 * 86400; % to convert volume transport from 10^6 m^3/s
                               % to m^3/yr
  for k=1:12,
    params.advect(k,:) = unitfac * (tr(k,:) + mix(k,:))/ params.volume(k);
  end
end

% initial PO4 distribution
po4_init = params.po4init;

% integrate
tspan = [0:50:3000];
po4 = ode23(@boxmodel_dgl_po4, tspan, po4_init);

% calculate RMS difference between final and initial PO4
po4_final = po4.y(:,end);
diff = po4_init - po4_final;
% f = sqrt( sum( (diff.^2).*params.volume ) ./ sum(params.volume) );
f = sqrt( sum( (diff.^2) ) );
%f = sum( abs(diff ) );

% prevent taking too strong mixing
if (max(abs(pvec))>40), 
    f = f + max(abs(pvec))-40;
end
if (min(pvec)<0), 
    f = f + 1000;
end

if (nargin>1),
    plot(po4.x,po4.y);
end

return
end
