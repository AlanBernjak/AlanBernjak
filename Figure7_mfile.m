% Script reproduces Figure7
% Corresponding data are in files 'Figure7_data

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% All AP time series are saved in the 'FabbriOutput' structure.

% The variable parameters are x_g_Kr (gKr), Ko and IHerg_shift_pa (IVshift)

clear
close all

x1 = load('Figure7_data');

allko=[];
allgkr=[];
allshift=[];
allCL=[];

BLCLthresh = 0.92;

tickgkr = [0.7 0.8 0.9 1];
tickKo = [3 3.5 4];

thresh = BLCLthresh;
% thresh = 0;

for st = 1:length(x1.parameters)
    
    allko = [allko; x1.parameters(st).Ko];
    allgkr = [allgkr; x1.parameters(st).x_g_Kr];
    allshift = [allshift; x1.parameters(st).IHerg_shift_pa];
    allCL = [allCL; x1.biomarkers(st).CL];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allCL = allCL-BLCLthresh;

ff1 = find(allCL>thresh);
ff2 = find(allCL<=thresh);

p11=plot3(allgkr(ff1),allshift(ff1),allko(ff1),'ro');
hold on
p12=plot3(allgkr(ff2),allshift(ff2),allko(ff2),'bo');    
    
set(p11,'markersize',4,'markeredgecolor','r','markerfacecolor','r');
set(p12,'markersize',4,'markeredgecolor','b','markerfacecolor','b');
axis tight

% view(-38,50)
view(126,25)

set(gca,'xlim',[0.69 1.01])
set(gca,'ylim',[-3.05 0.05])
set(gca,'zlim',[2.9 4.1])

set(gca,'xtick',tickgkr,'xticklabel',tickgkr)
set(gca,'ztick',tickKo,'zticklabel',tickKo)

xlabel('{\itg}_{Kr}')
yl=ylabel('IVshift (mV)');
% set(yl,'interpreter','tex')

zlabel('{\itK}o (mmol/L)')


ll=legend('CL>0.92s','CL<0.92s');

set(ll,'fontsize',12,'position',[0.70 0.83 0.21 0.11])

set(gca,'fontsize',14)

clear st

hold on
fill3([1 1 0.7 0.7],[0 0 -3 -3],[3 4 4 3],[0.8 0.8 0.8],'HandleVisibility','off')

% print -dtiff -r300dpi Ko3d_v3.tiff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

col1 = [0.8 0.8 0.8];
p1 = plot(allko,allCL,'o');

% set(p1,'markersize',4,'markeredgecolor',col1,'markerfacecolor',col1,'HandleVisibility','off')
set(p1,'markersize',4,'markeredgecolor',col1,'markerfacecolor',col1)
set(gca,'xlim',[2.95 4.05],'xtick',[3 3.2 3.4 3.6 3.8 4],'xticklabel',[3 3.2 3.4 3.6 3.8 4])
set(gca,'ylim',[0.74 1.35],'ytick',[0.8 0.9 1 1.1 1.2 1.3],'yticklabel',[0.8 0.9 1 1.1 1.2 1.3])

xlabel('{\itK}o (mmol/L)')
yl = ylabel(['CL (s)']);

%%%
setKo = 3:0.1:4;
setCL_hypo = [1.0362 1.0282 1.0200 1.0102 0.9991 0.9892 0.9771 0.9653 0.9535 0.9412 0.9285];
setCL_normg = [0.9392 0.9417 0.9433 0.9437 0.9432 0.9414 0.9386 0.9356 0.9318 0.9266 0.9209];

hold on
p3=plot(setKo,setCL_normg,'bo');
set(p3,'markersize',6,'markeredgecolor','k','markerfacecolor','b')

p2=plot(setKo,setCL_hypo,'ro');
set(p2,'markersize',6,'markeredgecolor','k','markerfacecolor','r')

%%%

set(gca,'fontsize',14)

ax = axis;
l1 = line([ax(1) ax(2)],[0.92 0.92]);
set(l1,'color','r','linewidth',2,'linestyle','--')

ll2=legend('{\itg}_{Kr}: [0.7 1], IVshift: [0 -3], {\itK}o: [3 4]','NG: {\itg}_{Kr}=1 & IVshift=0mV','LG: {\itg}_{Kr}=0.7 & IVshift=-3mV','Baseline ({\itK}o=4mmol/L): CL = 0.92s')
set(ll2,'position',[0.39 0.72 0.56 0.26],'fontsize',12)

% print -dtiff -r300dpi CLKo2d3.tiff
