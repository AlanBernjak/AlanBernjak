% Script reproduces Figure5
% Corresponding data are in files 'Figure5_data'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% For Figure5, parameters IHerg_shift_pa (IVshift) and x_g_Kr (gKr) 
% and biomarker CL have variable values, other parameters are constant 

clear
close all

x1 = load('Figure5_data');

allgkr=[];
allshift=[];
allCL=[];

allAPD = [];
allAPA = [];
allMDP1 = [];
allDDR = [];
allDDR100 = [];


% BLCLthresh = 0.69;
BLCLthresh = 0.92;
% vagCLthresh = 1.27;
% sympCLthresh = 0.69;


tickgkr = [0.7 0.75 0.8 0.85 0.9 0.95 1];

thresh = BLCLthresh;
% thresh = 0;

for st = 1:length(x1.parameters)
    
    allgkr = [allgkr; x1.parameters(st).x_g_Kr];
    allshift = [allshift; x1.parameters(st).IHerg_shift_pa];
    allCL = [allCL; x1.biomarkers(st).CL];
    
    allAPD = [allAPD; x1.biomarkers(st).APD90];
    allMDP1 = [allMDP1; x1.biomarkers(st).MDP1];
    allDDR = [allDDR; x1.biomarkers(st).DDR];
    allDDR100 = [allDDR100; x1.biomarkers(st).DDR100];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allCL = allCL-BLCLthresh;

ff1 = find(allCL>thresh);
ff2 = find(allCL<=thresh);

p11=plot(allgkr(ff1),allshift(ff1),'ro');
hold on
p12=plot(allgkr(ff2),allshift(ff2),'bo');    
    
set(p11,'markersize',3,'markerfacecolor','r');
set(p12,'markersize',3,'markerfacecolor','b');
axis tight

set(gca,'xlim',[0.69 1.01])
set(gca,'ylim',[-3.05 0.05])
set(gca,'xtick',tickgkr,'xticklabel',tickgkr)
xlabel('gKr')

set(gca,'fontsize',12)
view(180,90)

hold on
ll=line([1 0.7],[0 -3]);
set(ll,'color',[0.45 0.45 0.45],'linestyle','--','linewidth',2)
% print -dtiff -r300dpi 2dgKrV12_12.tiff

clear st

%%%%%%%%

f1=figure;

ff1_2 = (allCL>thresh);
ff2_2 = (allCL<=thresh);

p21=plot3(allgkr(ff1_2),allshift(ff1_2),allCL(ff1_2),'ro');
hold on
p22=plot3(allgkr(ff2_2),allshift(ff2_2),allCL(ff2_2),'bo');
view(155,31)
ll=legend('CL>0.92s','CL<0.92s');
set(ll,'position',[0.64 0.79 0.23 0.13],'fontsize',14)

ax = axis;
ll2 = plot3([1 0.7],[0 -3],[0.92 0.92],'HandleVisibility','off');
set(ll2,'color',[0.45 0.45 0.45],'linestyle','--','linewidth',2)

set(p21,'markersize',3,'markerfacecolor','r');
set(p22,'markersize',3,'markerfacecolor','b');
axis tight

set(gca,'xlim',[0.69 1.01])
set(gca,'ylim',[-3.05 0.05])
set(gca,'xtick',tickgkr,'xticklabel',tickgkr)
xlabel('{\it g}_{Kr}')
ylabel('IVshift (mV)')
zlabel('CL (s)')
% title('Ko: 4.0 mmol/l')

v1=[min(allCL(ff1_2)); max(allCL(ff1_2)); mean(allCL(ff1_2)); median(allCL(ff1_2))];
v2=[min(allCL(ff2_2)); max(allCL(ff2_2)); mean(allCL(ff2_2)); median(allCL(ff2_2))];
    
[v1 v2]

set(gca,'fontsize',14)

grid on

% print -dtiff -r300dpi 2dgKrV12_22_3.tiff

%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
sc=scatter(allgkr,allCL,[],allshift,'filled');

% color = [ones(200,1),(200-1:-1:0)'/199,zeros(200,1)];
color = [zeros(200,1),zeros(200,1),(200-1:-1:0)'/199];
colormap(color);

cb=colorbar;
set(cb,'Limits',[-3.01 0.01],'fontsize',14);
set(cb,'ticks',[-3 -1.5 0])
set(cb,'position',[0.22 0.63 0.03 0.26],'axislocation','in')

set(gca,'xlim',[0.69 1.01])
set(gca,'ylim',[0.74 1.21])

xlabel('{\it g}_{Kr}')
ylabel('CL (s)')

set(gca,'fontsize',14)

set(gca,'xtick',[0.7 0.8 0.9 1],'ytick',[0.8 0.9 1 1.1 1.2])

tb = annotation('textbox');
tb.FontSize = 14;
tb.Position = [0.2,0.91,0.27,0.07];
tb.EdgeColor = 'none';
tb.String = {'IVshift (mV)'};

ll=line([0.7 1],[0.92 0.92]);
set(ll,'color','r','linestyle','--','linewidth',2)

% print -dtiff -r300dpi 2dgKrV12_5.tiff

