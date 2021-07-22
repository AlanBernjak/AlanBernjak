% Script reproduces Figure6a
% Corresponding data are in files 'Figure6a_data1' and 'Figure6a_data2'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% For Figure5, parameters Ko and x_g_Kr (gKr) are variable 
% the other parameters are constant 

% Data for normoglycaemia and low glucose are saved in separate structures

close all
clear

x1=load('Figure6a_data1');

allgKr1 = [];
allCL1 = [];
allMDP1 = [];
allKo1 = [];

for st = 1:length(x1.parameters)
    
    allgKr1 = [allgKr1; x1.parameters(st).x_g_Kr];
    allCL1 = [allCL1; x1.biomarkers(st).CL];
    allMDP1 = [allMDP1; x1.biomarkers(st).MDP2];
    allKo1 = [allKo1; x1.parameters(st).Ko];
    
end

clear x1

%%%

x2=load('Figure6a_data2');

allgKr2 = [];
allCL2 = [];
allMDP2 = [];
allKo2 = [];

for st = 1:length(x2.parameters)
    
    allgKr2 = [allgKr2; x2.parameters(st).x_g_Kr];
    allCL2 = [allCL2; x2.biomarkers(st).CL];
    allMDP2 = [allMDP2; x2.biomarkers(st).MDP2];
    allKo2 = [allKo2; x2.parameters(st).Ko];
    
end

clear x2

%%%

figure

p1=plot(allKo1(1:45),allCL1(1:45),'bo');
hold on
p2=plot(allKo2(1:45),allCL2(1:45),'ro');

set(p1,'markersize',4,'markerfacecolor','b');
set(p2,'markersize',4,'markerfacecolor','r');

xlabel('{\itK}o (mmol/L)')
ylabel('CL (s)')

set(gca,'xlim',[0.5 5.1])
set(gca,'ylim',[0.73 1.61])

set(gca,'ytick',[0.8 1 1.2 1.4 1.6],'yticklabel',[0.8 1 1.2 1.4 1.6])
% set(gca,'xtick',[2 3 4 5],'xticklabel',[2 3 4 5])

ll = legend('Normoglycaemia','Low glucose');
set(ll,'fontsize',14)

ax = axis;
ax=axis;
l1=line([3 3],[ax(3) ax(4)],'HandleVisibility','off');
l2=line([4 4],[ax(3) ax(4)],'HandleVisibility','off');

set([l1 l2],'linewidth',2)

set(gca,'fontsize',14)
set(gca,'linewidth',1)

set(gcf,'position',[680 300 349 420])

% print -dtiff -r300dpi CLKoLowG3.tiff
