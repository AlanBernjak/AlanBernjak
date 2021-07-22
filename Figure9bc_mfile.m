% Script reproduces Figure9b and Figure9c
% Corresponding data are in files 'Figure9bc_data1','Figure9bc_data2' and 'Figure9bc_data3',

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% The script also needs function 'boxutil3.m' to plot Figure9c

clear
close all

x1 = load('Figure9bc_data1');

allko=[];
allgkr=[];
allshift=[];
allCL=[];

BLCLthresh = 0.92;
SympCLthresh = 0.69;
VagalCLthresh = 1.27;

tickgkr = [0.7 0.8 0.9 1];
tickKo = [3 3.5 4];

thresh = 1.27;
% thresh = 0;

for st = 1:length(x1.parameters)
    
    allko = [allko; x1.parameters(st).Ko];
    allgkr = [allgkr; x1.parameters(st).x_g_Kr];
    allshift = [allshift; x1.parameters(st).IHerg_shift_pa];
    allCL = [allCL; x1.biomarkers(st).CL];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = load('Figure9bc_data2');
x2 = load('Figure9bc_data3');

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

% ff1 = find(allCL>thresh);
% ff2 = find(allCL<=thresh);
% 
% p11=plot3(allgkr(ff1),allshift(ff1),allko(ff1),'ro');
% hold on
% p12=plot3(allgkr(ff2),allshift(ff2),allko(ff2),'bo');    
%     
% set(p11,'markersize',4,'markeredgecolor','r','markerfacecolor','r');
% set(p12,'markersize',4,'markeredgecolor','b','markerfacecolor','b');
% axis tight
% 
% % view(-38,50)
% view(126,25)
% axis tight
% 
% % set(gca,'xlim',[0.695 1.005])
% % set(gca,'ylim',[-3.05 0.05])
% % set(gca,'zlim',[2.9 4.1])
% 
% set(gca,'xtick',tickgkr,'xticklabel',tickgkr)
% set(gca,'ztick',tickKo,'zticklabel',tickKo)
% 
% xlabel('{\it g}_{Kr}')
% yl=ylabel('IVshift (mV)');
% % set(yl,'interpreter','tex')
% 
% zlabel('{\it K}o (mmol/L)')
% 
% 
% ll=legend('CL>1.27s','CL<1.27s');
% 
% set(ll,'fontsize',12,'position',[0.70 0.83 0.21 0.11]);
% 
% set(gca,'fontsize',12)
% 
% clear st
% 
% hold on
% fill3([1 1 0.7 0.7],[0 0 -3 -3],[3 4 4 3],[0.8 0.8 0.8],'HandleVisibility','off')
% 
% % print -dtiff -r300dpi Ko3d_v1.tiff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
% set(gcf,'position',[414 300 700 471])
set(gcf,'position',[914 352 520 471])


hold on

col1 = [0.8 0.8 0.8];
p1 = plot(allko,allCL,'o');

% set(p1,'markersize',4,'markeredgecolor',col1,'markerfacecolor',col1,'HandleVisibility','off')
set(p1,'markersize',4,'markeredgecolor',col1,'markerfacecolor',col1)

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

set(gca,'ylim',[0.85 4])

ax = axis;
l1 = line([ax(1) ax(2)],[0.92 0.92]);
set(l1,'color',[0.45 0.45 0.45],'linewidth',2,'linestyle','--');

l2 = line([ax(1) ax(2)],[1.27 1.27]);
set(l2,'color','r','linewidth',2,'linestyle','--');

ll2=legend('{\itg}_{Kr}: [0.7 1], IVshift: [0 -3], {\itK}o: [3 4]',...
           'NG: {\itg}_{Kr}=1 & IVshift=0',...
           'LG: {\itg}_{Kr}=0.7 & IVshift=-3mV');
       
set(ll2,'position',[0.4 0.76 0.4 0.17],'fontsize',12);

ttext = text(3.28,4.18,'Parasymp.Activity +');
set(ttext,'fontsize',13)

set(gca,'linewidth',1)

set(gca,'fontsize',14)
% set(gcf,'position',[914 352 620 471])

% print -dtiff -r300dpi CLKoVag4.tiff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf,'position',[914 352 120 471])

% axes(sp2);
hold on
% boxplot(log10(allCL))
boxutil3_v2(log10(allCL),1,1,0.3,'o',1,1.5,[0.8 0.8 0.8])
ylabel('CL (s)')

yt = [0 0.18 0.3 0.48 0.7 1 1.4];
yt2 = [1 1.5 2 3 5 10 25];
box('off')

set(gca,'ytick',yt,'yticklabel',yt2)
set(gca,'xtick',[])
set(gca,'linewidth',1)

set(gca,'xlim',[0.5 1.5])
set(gca,'fontsize',14)
set(gca,'linewidth',1)

% print -dtiff -r300dpi CLKoVag42.tiff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















