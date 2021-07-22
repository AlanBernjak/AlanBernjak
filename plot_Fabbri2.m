% plots the CL data in the Ko,gKr,IVshift parameter space
% input parameters and AP features data are saved in the testFabbri2 structure

clear
close all

x1 = load('testFabbri2');

allko=[];
allgkr=[];
allshift=[];
allCL=[];

tickgkr = [0.7 0.8 0.9 1];
tickKo = [3 3.5 4];

for st = 1:length(x1.parameters)
    
    allko = [allko; x1.parameters(st).Ko];
    allgkr = [allgkr; x1.parameters(st).x_g_Kr];
    allshift = [allshift; x1.parameters(st).IHerg_shift_pa];
    allCL = [allCL; x1.biomarkers(st).CL];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

p11=plot3(allgkr,allshift,allko,'bo');
    
set(p11,'markersize',4,'markeredgecolor','r','markerfacecolor','r');
% axis tight

% view(-38,50)
view(126,25)
% axis tight

set(gca,'xtick',tickgkr,'xticklabel',tickgkr)
set(gca,'ztick',tickKo,'zticklabel',tickKo)

xlabel('{\it g}_{Kr}')
yl=ylabel('IVshift (mV)');
% set(yl,'interpreter','tex')

zlabel('{\it K}o (mmol/L)')

set(gca,'fontsize',12)

clear st

% print -dtiff -r300dpi Ko3d_v1.tiff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



