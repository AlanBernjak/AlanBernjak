% Script reproduces Figure12
% Corresponding data are in files 'Figure12_data1' and 'Figure12_data2'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% AP time series and currents are saved in the 'SeveriOutput' structure.
% All time series are cropped and only data displayed in the figure are
% available. The time vector shows the actual time of the simulation.


clear 
close all

load Figure12_data1

f1=figure;
% set(f1,'position',[178 87 668 900])

% ind1 = 700000;
% ind2 = ind1+700;

time = SeveriOutput.time;
volt = SeveriOutput.volt;

% time = (time-450)*100;

If = SeveriOutput.If;
INaCa = SeveriOutput.INaCa;
ICaL = SeveriOutput.ICaL;
ICaT = SeveriOutput.ICaT;
IKr = SeveriOutput.IKr;

%%%%%%%%%

% If = If(ind1:ind2);
% INaCa = INaCa(ind1:ind2);
% ICaL = ICaL(ind1:ind2);
% ICaT = ICaT(ind1:ind2);
% IKr = IKr(ind1:ind2);
% 
% biomarkers = Severi_findfiducial(time,volt);
% CL1 = biomarkers.CL;
% 
% time = time(ind1:ind2);
% volt = volt(ind1:ind2);

%%%%%%%%%

sp(1)=subplot(211);
p1=plot(time,volt);
% axis tight

% set(p1,'linewidth',2)



% t1=title(['Rabbit SAN, gKr=1, CL: 355 ms']);
yl1 = ylabel('Voltage (mV)');

set(gca,'ylim',[-65 -25])
set(gca,'xlim',[492.48 492.8])

tx = text(492.55,-30,'{\itg}_{Kr}=1.0, CL=0.36s');
set(tx,'fontsize',18);

set(gca,'fontsize',18)
set(gca,'linewidth',1)
set(gca,'xticklabel',[])

pos1=get(gca,'position');
set(gca,'position',[0.13 0.79 0.81 0.19]);
set(gca,'box','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp(2)=subplot(212);
hold on
p2=plot(time,If);
p3=plot(time,INaCa);
p4=plot(time,ICaL);
p5=plot(time,ICaT);
p6=plot(time,IKr);
% axis tight

yl2 = ylabel('Current (pA/pF)');
xl2 = xlabel('Time (s)');
set(sp(2),'ylim',[-0.3 0.2])

set([p1,p2,p3,p4,p5,p6],'linewidth',2)

set(gca,'position',[0.13 0.11 0.81 0.67])

leg1 = legend('I_f','I_{NaCa}','I_{CaL}','I_{CaT}','I_{Kr}');
set(leg1,'position',[0.16 0.53 0.13 0.168])

set(gca,'xlim',[492.48 492.8])

linkaxes(sp,'x')

set(gca,'fontsize',18)
set(gca,'linewidth',1)

set(gca,'ytick',[-0.3 -0.2 -0.1 0 0.1 0.2])


pos2=get(gca,'position');
set(gca,'position',[0.13 0.097 0.81 0.67]);

% print -dtiff SeveriCurr_gKr1x3.tiff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Figure12_data2
% IKr205 = FabbriOutput.IKr;

f2=figure;
% set(f2,'position',[977 87 668 900])

% ind1 = 653200;
% ind2 = ind1+700;

time2 = SeveriOutput.time;
volt2 = SeveriOutput.volt;

If2 = SeveriOutput.If;
INaCa2 = SeveriOutput.INaCa;
ICaL2 = SeveriOutput.ICaL;
ICaT2 = SeveriOutput.ICaT;
IKr2 = SeveriOutput.IKr;

%%%%%%%%%%%%%%%%%%

% If2 = If2(ind1:ind2);
% INaCa2 = INaCa2(ind1:ind2);
% ICaL2 = ICaL2(ind1:ind2);
% ICaT2 = ICaT2(ind1:ind2);
% IKr2 = IKr2(ind1:ind2);
% 
% biomarkers2 = Severi_findfiducial(time2,volt2);
% CL2 = biomarkers2.CL;
% 
% time2 = time2(ind1:ind2);
% volt2 = volt2(ind1:ind2);

%%%%%%%%%%%%%%%%%%%%%%%%%%5

sp2(1)=subplot(211);
p1=plot(time2,volt2);

yl3 = ylabel('Voltage (mV)');
% xl2 = xlabel('Time (ms)')

% axis tight



% t2=title(['Rabbit SAN, gKr=0.5, CL: 411 ms']);

set(gca,'ylim',[-65 -25])

tx = text(498.47,-30,'{\itg}_{Kr}=0.5, CL=0.41s');
set(tx,'fontsize',18);

set(gca,'fontsize',18)
set(gca,'linewidth',1)
set(gca,'xticklabel',[])

pos1=get(gca,'position');
set(gca,'position',[0.13 0.79 0.81 0.19]);
set(gca,'box','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sp2(2)=subplot(212);
hold on
p2=plot(time2,If2);
p3=plot(time2,INaCa2);
p4=plot(time2,ICaL2);
p5=plot(time2,ICaT2);
p6=plot(time2,IKr2);
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

set(gca,'ytick',[-0.3 -0.2 -0.1 0 0.1 0.2])

set(gca,'xlim',[498.38 498.68])

% print -dtiff SeveriCurr_gKr05x3.tiff
