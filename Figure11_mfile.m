% Function reproduces Figure11a and Figure11b
% Corresponding data are in files 'Figure11_data1' to 'Figure11_data6'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% All AP time series are saved in the 'FabbriOutput' structure.


function Figure11_mfile

clear
close all

load Figure11_data1

time1 = SeveriOutput.time;
volt1 = SeveriOutput.volt;
CL1 = biomarkers.CL;

load Figure11_data2

time2 = SeveriOutput.time;
volt2 = SeveriOutput.volt;
CL2 = biomarkers.CL;

load Figure11_data3

time3 = SeveriOutput.time;
volt3 = SeveriOutput.volt;
CL3 = biomarkers.CL;

load Figure11_data4

time4 = SeveriOutput.time;
volt4 = SeveriOutput.volt;
CL4 = biomarkers.CL;

load Figure11_data5

time5 = SeveriOutput.time;
volt5 = SeveriOutput.volt;
CL5 = biomarkers.CL;

load Figure11_data6

time6 = SeveriOutput.time;
volt6 = SeveriOutput.volt;
CL6 = biomarkers.CL;

chosenp = 3;

dpos1 = selectsection(time1,volt1,chosenp);
reltime1=time1-time1(dpos1);

f1=figure;
p1=plot(reltime1(dpos1:end),volt1(dpos1:end));

dpos2 = selectsection(time2,volt2,chosenp);
reltime2=time2-time2(dpos2);
hold on
p2=plot(reltime2(dpos2:end),volt2(dpos2:end));

dpos3 = selectsection(time3,volt3,chosenp);
reltime3=time3-time3(dpos3);
hold on
p3=plot(reltime3(dpos3:end),volt3(dpos3:end));

dpos4 = selectsection(time4,volt4,chosenp);
reltime4=time4-time4(dpos4);
hold on
p4=plot(reltime4(dpos4:end),volt4(dpos4:end));

dpos5 = selectsection(time5,volt5,chosenp);
reltime5=time5-time5(dpos5);
hold on
p5=plot(reltime5(dpos5:end),volt5(dpos5:end));

posx = 400000;
dpos6 = selectsection(time6,volt6,chosenp);
reltime6=time6-time6(posx);
hold on
p6=plot(reltime6(posx:end),volt6(posx:end));

set(p1,'linewidth',3)
set([p2,p3,p4,p5,p6],'linewidth',1.5)

set(gca,'xlim',[-0.01 1.01])
set(gca,'linewidth',1)

xl=xlabel('Time (s)');
yl=ylabel('Voltage (mV)');

legend({'{\it g}_{Kr} = 1.0','{\it g}_{Kr} = 0.8','{\it g}_{Kr} = 0.6','{\it g}_{Kr} = 0.4','{\it g}_{Kr} = 0.2','{\it g}_{Kr} = 0'})

set(gca,'fontsize',14)
set(gca,'box','off')

allgkr = [1 0.8 0.6 0.4 0.2];
allCL = [CL1 CL2 CL3 CL4 CL5];


% print -dtiff -r300dpi APgKr_Severi2.tiff

%%%%%%%%%%%%%%%%

f2=figure;
p1=plot(allgkr,allCL,'bo');
set(p1,'markersize',6,'markerfacecolor','b');


set(gca,'xlim',[0.16 1.04])
xl=xlabel('{\it g}_{Kr}');
yl=ylabel('CL (s)');

set(gca,'fontsize',14,'box','off')

% set(gcf,'position',[680 558 351 420])

set(gca,'linewidth',1)

% print -dtiff -r300dpi CLgKr_Severi3.tiff


function dpos = selectsection(time1,volt1,chosenp)

[peaks1, peakpos1]=findpeaks(volt1);

chosenp = 3;

beatn = length(peaks1)-chosenp;

dvolt1 = diff(volt1);

[dpeaks1, dpeakpos1]=findpeaks(dvolt1,'MinPeakDistance',300);

dpos = dpeakpos1(end-chosenp-1);

reltime1=time1-time1(dpos);
