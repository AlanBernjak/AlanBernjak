% Script reproduces Figure3a
% Corresponding data are in files 'Figure3a_data1' and 'Figure3a_data2'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% For Figure3a, parameter x_g_Kr and biomarker CL have 
% variable values, other parameters are constant 

close all
clear

x1=load('Figure3a_data1');

allgKr = [];
allCL = [];
allMDP = [];

for st = 1:length(x1.parameters)
    
    allgKr = [allgKr; x1.parameters(st).x_g_Kr];
    allCL = [allCL; x1.biomarkers(st).CL];
    allMDP = [allMDP; x1.biomarkers(st).MDP2];
    
end

x2=load('Figure3a_data2');

for st = 5%1:length(x2.parameters)
    
    allgKr = [allgKr; x2.parameters(st).x_g_Kr];
    allCL = [allCL; x2.biomarkers(st).CL];
    allMDP = [allMDP; x2.biomarkers(st).MDP2];
    
end

figure

p1=plot(allgKr,allCL,'bo');
hold on

set(p1,'markersize',4,'markerfacecolor','b');

xl = xlabel('{\it g}_{Kr}');
ylabel('CL (s)')

set(gca,'xlim',[0.07 1.02])
set(gca,'ylim',[0.44 0.95])

set(gca,'ytick',[0.5 0.6 0.7 0.8 0.9],'yticklabel',[0.5 0.6 0.7 0.8 0.9])

ax = axis;
ax=axis;
l1=line([0.7 0.7],[ax(3) ax(4)]);
l2=line([1 1],[ax(3) ax(4)]);

set([l1 l2],'linewidth',2)

set(gca,'fontsize',14)
set(gca,'linewidth',1)

% set(gcf,'position',[680 558 351 420])

% print -dtiff CLgKr2.tiff
