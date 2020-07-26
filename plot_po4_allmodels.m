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

figure(2)
clf
h0 = plot(po4_woa, po4_cl,'.');
set(h0,'LineWidth',2,'MarkerSize',50,'Color',[0.8 0.2 0.1]);
hold on 
h1 = plot(po4_woa, po4_1l,'.');
set(h1,'LineWidth',2,'MarkerSize',40,'Color',[0 0.15 0.70]);
h2 = plot(po4_woa, po4_2l1,'.');
%set(h2,'LineWidth',2,'MarkerSize',30,'Color',[0.3660 0.5740 0.0880]);
%set(h2,'LineWidth',2,'MarkerSize',30,'Color',[0.16 0.80 0.05]);
set(h2,'LineWidth',2,'MarkerSize',30,'Color',[0.48 0.75 0.26]);
h3 = plot(po4_woa, po4_2l2,'.');
% set(h3,'LineWidth',2,'MarkerSize',20,'Color',[0.8290 0.5940 0.0250]);
set(h3,'LineWidth',2,'MarkerSize',20,'Color',[0.960 0.80 0.37]);
%set(h3,'LineWidth',2,'MarkerSize',20,'Color',[0.8290 0.040 0.850]);
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

print('PO4_all4boxmodels.png','-dpng','-r600');
