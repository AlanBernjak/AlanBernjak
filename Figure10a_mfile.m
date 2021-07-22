% Script reproduces Figure10a
% Corresponding data are in files 'Figure10a_data1','Figure10a_data2' and 'Figure10a_data3',

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 


clear
close all

x3 = load('Figure10a_data1');

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

for st = 1:length(x3.parameters)
    
    allko = [allko; x3.parameters(st).Ko];
    allgkr = [allgkr; x3.parameters(st).x_g_Kr];
    allshift = [allshift; x3.parameters(st).IHerg_shift_pa];
    allCL = [allCL; x3.biomarkers(st).CL];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = load('Figure10a_data2');
x2 = load('Figure10a_data3');

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

ff1 = find((allCL>thresh)&(allCL<2));
ff2 = find((allCL<=thresh)&(allCL>0.5));
ff3 = find(allCL>=2);
ff4 = find(allCL<=0.5);

allcol = get(0,'defaultAxesColorOrder');

p11=plot3(allgkr(ff1),allshift(ff1),allko(ff1),'ro'); 
hold on
p12=plot3(allgkr(ff2),allshift(ff2),allko(ff2),'bo');
p13=plot3(allgkr(ff3),allshift(ff3),allko(ff3),'go');   
p14=plot3(allgkr(ff4),allshift(ff4),allko(ff4),'o');   
    
set(p11,'markersize',3,'markerfacecolor','r');
set(p12,'markersize',3,'markerfacecolor','b');
set(p13,'markersize',3,'markerfacecolor','g');
set(p14,'markersize',3,'markerfacecolor',allcol(4,:),'markeredgecolor',allcol(4,:));
axis tight

% view(-38,50)
view(127,29)

set(gca,'xlim',[0.69 1.01])
set(gca,'ylim',[-3.05 0.05])
set(gca,'xtick',tickgkr,'xticklabel',tickgkr)
xlabel('{\itg}_{Kr}')
ylabel('IVshift')
zlabel('{\itK}o (mmol/L)')

ll=legend('CL > 1.27s','CL < 1.27s','CL > 2s','CL = N/A')

clear st

set(ll,'fontsize',12,'position',[0.73 0.72 0.23 0.21]);

set(gca,'fontsize',14)

clear st

hold on
fill3([1 1 0.7 0.7],[0 0 -3 -3],[3 4 4 3],[0.8 0.8 0.8],'HandleVisibility','off')

