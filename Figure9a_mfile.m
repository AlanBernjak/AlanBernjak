% Script reproduces Figure9a
% Corresponding data are in files 'Figure9a_data1','Figure9a_data2' and 'Figure9a_data3',

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

clear
close all

x1 = load('Figure9a_data1');

allko=[];
allgkr=[];
allshift=[];
allCL=[];

BLCLthresh = 0.92;
SympCLthresh = 0.69;
VagalCLthresh = 1.27;

tickgkr = [0.7 0.8 0.9 1];
tickKo = [3 3.5 4];

thresh = 0.69;
% thresh = 0;

for st = 1:length(x1.parameters)
    
    allko = [allko; x1.parameters(st).Ko];
    allgkr = [allgkr; x1.parameters(st).x_g_Kr];
    allshift = [allshift; x1.parameters(st).IHerg_shift_pa];
    allCL = [allCL; x1.biomarkers(st).CL];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = load('Figure9a_data2');
x2 = load('Figure9a_data3');

allKo1=[];
allKo2=[];
allCL1=[];
allCL2=[];

for st2 = 1:length(x1.parameters)        
    allCL1 = [allCL1; x1.biomarkers(st2).CL];    
    allCL2 = [allCL2; x2.biomarkers(st2).CL];    
    allKo1 = [allKo1; x1.parameters(st2).Ko];    
    allKo2 = [allKo2; x2.parameters(st2).Ko];    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
axis tight

set(gca,'xlim',[0.695 1.005])
set(gca,'ylim',[-3.05 0.05])
set(gca,'zlim',[2.9 4.1])

set(gca,'xtick',tickgkr,'xticklabel',tickgkr)
set(gca,'ztick',tickKo,'zticklabel',tickKo)

xlabel('{\it g}_{Kr}')
yl=ylabel('IVshift (mV)');
% set(yl,'interpreter','tex')

zlabel('{\it K}o (mmol/L)')


ll=legend('CL>0.69s','CL<0.69s');

set(ll,'fontsize',12,'position',[0.70 0.83 0.21 0.11]);

set(gca,'fontsize',12)

clear st

hold on
fill3([1 1 0.7 0.7],[0 0 -3 -3],[3 4 4 3],[0.8 0.8 0.8],'HandleVisibility','off')

% print -dtiff -r300dpi Ko3d_v1.tiff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

col1 = [0.8 0.8 0.8];
p1 = plot(allko,allCL,'o');

set(p1,'markersize',4,'markeredgecolor',col1,'markerfacecolor',col1)
% set(p1,'markersize',4,'markeredgecolor',col1,'markerfacecolor',col1,'HandleVisibility','off')
set(gca,'xlim',[2.95 4.05],'xtick',[3 3.2 3.4 3.6 3.8 4],'xticklabel',[3 3.2 3.4 3.6 3.8 4])
% set(gca,'ylim',[0.75 1.2533],'ytick',[0.8 0.9 1 1.1 1.2 1.3],'yticklabel',[0.8 0.9 1 1.1 1.2 1.3])

xlabel('{\it K}o (mmol/L)')
yl = ylabel(['CL (s)']);

%%%

hold on
p3=plot(allKo1,allCL1,'bo');
set(p3,'markersize',6,'markeredgecolor','k','markerfacecolor','b')

p2=plot(allKo1,allCL2,'ro');
set(p2,'markersize',6,'markeredgecolor','k','markerfacecolor','r')

%%%

% set(gca,'fontsize',12,'xdir','reverse')

% set(gca,'ylim',[0.45 0.96])
set(gca,'ylim',[0.51 1.1])

ax = axis;
l1 = line([ax(1) ax(2)],[0.92 0.92]);
set(l1,'color',[0.45 0.45 0.45],'linewidth',2,'linestyle','--');

l2 = line([ax(1) ax(2)],[0.69 0.69]);
set(l2,'color','r','linewidth',2,'linestyle','--');

% ll2=legend('{Symp.Act.+\itg}_{Kr}: [0.7 1], IVshift: [0 -3], {\itK}o: [3 4]',...
%            'NG+Symp.Act.: {\itg}_{Kr}=1 & IVshift=0',...
%            'LG+Symp.Act.: {\itg}_{Kr}=0.7 & IVshift=-3mV',...
%            'Baseline ({\itK}o=4mmol/L): CL=0.92s',...
%            'Baseline + Symp.Act.: CL = 0.69s');
ll2=legend('{\itg}_{Kr}: [0.7 1], IVshift: [0 -3], {\itK}o: [3 4]',...
           'NG: {\itg}_{Kr}=1 & IVshift=0',...
           'LG: {\itg}_{Kr}=0.7 & IVshift=-3mV');    
       
set(ll2,'position',[0.38 0.73 0.47 0.21],'fontsize',12);

ttext = text(3.3,1.14,'Symp.Activity +');
set(ttext,'fontsize',13)

set(gca,'fontsize',14)
set(gca,'linewidth',1)

set(gcf,'position',[914 352 520 471])

box off
% print -dtiff -r300dpi CLKoSymp4.tiff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



















