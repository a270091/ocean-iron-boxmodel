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

% extract PO4 from the model runs
po4_cl  = a0(:,1);
po4_1l  = a(:,1);
po4_2l1 = b(:,1);
po4_2l2 = c(:,1);
po4_3l1 = d(:,1);
po4_3l2 = d2(:,1);

% read box model initial parameters, to have the box names and the 
% PO4 distribution from WOA
global params
boxmodel_init_params
po4_woa = params.po4init;

figure(2)
clf
h0 = plot(po4_woa, po4_cl,'.');
set(h0,'LineWidth',2,'MarkerSize',60,'Color',[0.8 0.2 0.1]);
hold on 
h1 = plot(po4_woa, po4_1l,'.');
set(h1,'LineWidth',2,'MarkerSize',50,'Color',[0 0.15 0.70]);
h4 = plot(po4_woa, po4_3l1,'.');
%set(h4,'MarkerSize',40, 'LineWidth',2,'Color',[0.00 0.75 0.00]);
set(h4,'MarkerSize',40, 'LineWidth',2,'Color',[120/256 221/256 120/256]);
h5 = plot(po4_woa, po4_3l2,'.');
set(h5,'MarkerSize',30, 'LineWidth',2,'Color',[62/256 150/256 55/256]);
% set(h5,'MarkerSize',30, 'LineWidth',2,'Color',[0.70 0.0 0.95]);
% h2 = plot(po4_woa, po4_2l1,'.');
% set(h2,'MarkerSize',30, 'LineWidth',2,'Color',[0.0 0.9 0.9]);
% set(h2,'LineWidth',2,'MarkerSize',30,'Color',[0.48 0.75 0.26]);
h3 = plot(po4_woa, po4_2l2,'.');
% set(h3,'LineWidth',2,'MarkerSize',20,'Color',[0.960 0.80 0.37]);
set(h3,'MarkerSize',20, 'LineWidth',2,'Color',[104/256 225/256 238/256]);
% set(h3,'MarkerSize',20, 'LineWidth',2,'Color',[0.90 0.90 0.0]);
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',12);
xlabel('WOA-derived PO_4');
ylabel('model PO_4');
yl = get(gca,'YLim');
hd = plot(yl,yl,'k');

% hand-edited annotations
for k=1:12
    xpos = [po4_woa(k),po4_woa(k)];
    if (k==1)
        yp1 = min([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.5*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==2)
        yp1 = min([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.7*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==3)
        yp1 = min([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.6*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==4)
        yp1 = min([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.4*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==5)
        yp1 = min([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.2*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==6)
        yp1 = max([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.3*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==7)
        yp1 = max([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.2*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==8)
        yp1 = max([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.8*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==9)
        yp1 = max([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.4*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==10)
        yp1 = min([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.8*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==11)
        yp1 = max([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.6*yl(2);
        ypos = ([yp1,yp2]);
    elseif (k==12)
        yp1 = min([po4_cl(k),po4_1l(k),po4_2l1(k),po4_2l2(k)]);
        yp2 = 0.9*yl(2);
        ypos = ([yp1,yp2]);
    else
        ypos = yl;
    end
    ha = plot(xpos,ypos,'k');
    ht = text(xpos(2),ypos(2),params.names{k});
    if (ypos(2)>xpos(2))
        set(ht,'FontSize',12 ,'Margin',5,'horizontalAlignment','right',...
            'verticalAlignment','bottom')
    else
        set(ht,'FontSize',12 ,'Margin',5,'horizontalAlignment','left',...
            'verticalAlignment','top')
    end
end

% add a legend
x0 = 2.6;
y0 = 0.85;
dy = 0.15;
dx = 0.15;
hl0 = plot(x0,y0,'.');
ht0 = text(x0+dx,y0,'CL');
set(hl0,'MarkerSize',60, 'LineWidth',2,'Color',[0.8 0.2 0.1]);
set(ht0,'FontSize',12)
hl1 = plot(x0,y0-dy,'.');
ht1 = text(x0+dx,y0-dy,'1L');
set(hl1,'MarkerSize',50, 'LineWidth',2,'Color',[0 0.15 0.70]);
set(ht1,'FontSize',12)
hl5 = plot(x0,y0-2*dy,'.');
ht5 = text(x0+dx,y0-2*dy,'1H');
%set(hl5,'MarkerSize',40, 'LineWidth',2,'Color',[0.70 0.0 0.95]);
set(hl5,'MarkerSize',40, 'LineWidth',2,'Color',[62/256 150/256 55/256]);
set(ht5,'FontSize',12)
hl4 = plot(x0,y0-3*dy,'.');
ht4 = text(x0+dx,y0-3*dy,'2H');
set(hl4,'MarkerSize',30, 'LineWidth',2,'Color',[120/256 221/256 120/256]);
set(ht4,'FontSize',12)
hl3 = plot(x0,y0-4*dy,'.');
ht3 = text(x0+dx,y0-4*dy,'2L');
% set(hl3,'MarkerSize',20, 'LineWidth',2,'Color',[0.90 0.90 0.00]);
set(hl3,'MarkerSize',20, 'LineWidth',2,'Color',[104/256 225/256 238/256]);
set(ht3,'FontSize',12)


print('PO4_allmodels_mar21.png','-dpng','-r600');

