% Script reproduces Figure4a
% Corresponding data are in files 'Figure4a_data'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% For Figure3a, parameter IHerg_shift_pa and biomarker CL have 
% variable values, other parameters are constant 

close all
clear

x1=load('Figure4a_data');

% x1=load('test20');

allgKr = [];
allCL = [];
allMDP = [];
allV12 = [];

for st = 1:length(x1.parameters)
    
    allgKr = [allgKr; x1.parameters(st).x_g_Kr];
    allCL = [allCL; x1.biomarkers(st).CL];
    allMDP = [allMDP; x1.biomarkers(st).MDP2];
    allV12 = [allV12; x1.parameters(st).IHerg_shift_pa];
    
end

figure

p1=plot(allV12,allCL,'bo');
hold on

set(p1,'markersize',4,'markerfacecolor','b');

xlabel('IVshift (mV)')
ylabel('CL (s)')

set(gca,'xlim',[-5.4 0.2])
set(gca,'ylim',[0.8 2.6])

set(gca,'ytick',[1 1.5 2 2.5 3],'yticklabel',[1 1.5 2 2.5 3])
set(gca,'xtick',[-5 -4 -3 -2 -1 0],'xticklabel',[-5 -4 -3 -2 -1 0])

ax = axis;
ax=axis;
l1=line([0 0],[ax(3) ax(4)]);
l2=line([-3 -3],[ax(3) ax(4)]);

set([l1 l2],'linewidth',2)

set(gca,'fontsize',14)
set(gca,'linewidth',1)

% set(gcf,'position',[680 558 351 420])

% print -dtiff CLV12_3.tiff
