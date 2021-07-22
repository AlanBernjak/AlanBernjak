% Script reproduces Figure13
% Corresponding data are in files 'Figure13_data1' and 'Figure13_data2'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% All AP time series and currents are saved in the 'FAbbriOutput' structure.

% The script also needs function 'Severi_findfiducial.m' to plot the Figures

clear 
close all

load('Figure13_data1')

f1=figure;
% set(f1,'position',[178 87 668 900])

ind1 = 234000;
ind2 = ind1+1300;

time = FabbriOutput.time;
volt = FabbriOutput.volt;

If = FabbriOutput.If;
INaCa = FabbriOutput.INaCa;
ICaL = FabbriOutput.ICaL;
ICaT = FabbriOutput.ICaT;
IKr = FabbriOutput.IKr;

sp(1)=subplot(211);
p1=plot(time(ind1:ind2),volt(ind1:ind2));
% axis tight

biomarkers = Severi_findfiducial(time,volt);
CL1 = biomarkers.CL;

% t1=title(['Human SAN, gKr=1, CL: 921 ms']);
yl1 = ylabel('Voltage (mV)');

set(gca,'ylim',[-65 -25])
set(gca,'xlim',[198.5 199.43])

tx = text(198.7,-30,'{\itg}_{Kr}=1, CL=0.92s');
set(tx,'fontsize',18);

set(gca,'fontsize',18)
set(gca,'linewidth',1)
set(gca,'xticklabel',[])

pos1=get(gca,'position');
set(gca,'position',[0.13 0.79 0.81 0.19]);
set(gca,'box','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp(2)=subplot(212);
hold on
p2=plot(time(ind1:ind2),If(ind1:ind2));
p3=plot(time(ind1:ind2),INaCa(ind1:ind2));
p4=plot(time(ind1:ind2),ICaL(ind1:ind2));
p5=plot(time(ind1:ind2),ICaT(ind1:ind2));
p6=plot(time(ind1:ind2),IKr(ind1:ind2));
% axis tight

yl2 = ylabel('Current (pA/pF)');
xl2 = xlabel('Time (s)');
set(sp(2),'ylim',[-0.3 0.2])

set([p1,p2,p3,p4,p5,p6],'linewidth',2)

set(gca,'position',[0.13 0.11 0.81 0.67])

leg1 = legend('I_f','I_{NaCa}','I_{CaL}','I_{CaT}','I_{Kr}');
set(leg1,'position',[0.16 0.53 0.11 0.168])

set(gca,'xlim',[198.5 199.43])
set(gca,'ytick',[-0.3 -0.2 -0.1 0 0.1 0.2])

linkaxes(sp,'x')

set(gca,'fontsize',18)
set(gca,'linewidth',1)

pos2=get(gca,'position');
set(gca,'position',[0.13 0.097 0.81 0.67]);

% print -dtiff FabbriCurr_gKr1x3.tiff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

load('Figure13_data2')
% IKr205 = FabbriOutput.IKr;

f2=figure;
% set(f2,'position',[977 87 668 900])

ind1 = 234000;
ind2 = ind1+1300;

time2 = FabbriOutput.time;
volt2 = FabbriOutput.volt;

If2 = FabbriOutput.If;
INaCa2 = FabbriOutput.INaCa;
ICaL2 = FabbriOutput.ICaL;
ICaT2 = FabbriOutput.ICaT;
IKr2 = FabbriOutput.IKr;

sp2(1)=subplot(211);
p1=plot(time2(ind1:ind2),volt2(ind1:ind2));
% axis tight

biomarkers2 = Severi_findfiducial(time2,volt2);
CL2 = biomarkers2.CL;

yl1 = ylabel('Voltage (mV)');

set(gca,'ylim',[-65 -25])
set(gca,'xlim',[192.6 193.43])

tx = text(192.75,-30,'{\itg}_{Kr}=0.5, CL=0.60s');
set(tx,'fontsize',18);

set(gca,'fontsize',18)
set(gca,'linewidth',1)
set(gca,'xticklabel',[])

pos1=get(gca,'position');
set(gca,'position',[0.13 0.79 0.81 0.19]);
set(gca,'box','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp2(2)=subplot(212);
hold on
p2=plot(time2(ind1:ind2),If2(ind1:ind2));
p3=plot(time2(ind1:ind2),INaCa2(ind1:ind2));
p4=plot(time2(ind1:ind2),ICaL2(ind1:ind2));
p5=plot(time2(ind1:ind2),ICaT2(ind1:ind2));
p6=plot(time2(ind1:ind2),IKr2(ind1:ind2));
% plot(time2(ind1:ind2),IKr205(ind1:ind2))
% axis tight

yl4 = ylabel('Current (pA/pF)');
xl2 = xlabel('Time (s)')

set([p1,p2,p3,p4,p5,p6],'linewidth',2)

set(sp2(2),'ylim',[-0.3 0.2])

set(gca,'fontsize',18)
set(gca,'linewidth',1)

pos1=get(gca,'position');
set(gca,'position',[0.13 0.097 0.81 0.67]);

set(gca,'ytick',[-0.3 -0.2 -0.1 0 0.1 0.2])

linkaxes(sp2,'x')

leg2 = legend('I_f','I_{NaCa}','I_{CaL}','I_{CaT}','I_{Kr}');

set(leg2,'position',[0.16 0.53 0.13 0.168])

set(gca,'xlim',[192.6 193.43])

% print -dtiff FabbriCurr_gKr05x3.tiff
