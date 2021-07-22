% Script reproduces Figure10b
% Corresponding data are in files 'Figure10b_data1' to 'Figure10b_data7'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% All AP time series are saved in the 'FabbriOutput' structure.


clear
close all

d1=load('Figure10b_data1');
d2=load('Figure10b_data2');
d3=load('Figure10b_data3');
d4=load('Figure10b_data4');
d5=load('Figure10b_data5');
d6=load('Figure10b_data6');
d7=load('Figure10b_data7');
% d8=load('test20x4');

x0 = 0.11;
y0 = 0.13;

dx = 0.85;
nfig = 7;
gapy = 0.02;

dy=(0.98-y0-nfig*gapy)/nfig;

ydim1=-85;
ydim2=40;

xdim1=0;
xdim2=42;

xt = [0 10 20 30 40];
yt = [-60 30];
fsize = 14;
lwid = 1.5;

figure
sp1 = subplot(711);
set(sp1,'position',[x0 y0+(nfig-1)*(gapy+dy) dx dy])
p1=plot(d1.FabbriOutput.time,d1.FabbriOutput.volt);
set(p1,'linewidth',lwid)
set(gca,'xticklabel',[],'ylim',[ydim1 ydim2],'xlim',[xdim1 xdim2],'xtick',xt,'ytick',yt)
box('off')
set(gca,'fontsize',fsize)

sp2 = subplot(712);
set(sp2,'position',[x0 y0+(nfig-2)*(gapy+dy) dx dy])
[i j]=size(d2.FabbriOutput);
p2=plot(d2.FabbriOutput(j).time,d2.FabbriOutput(j).volt);
set(p2,'linewidth',lwid)
set(gca,'xticklabel',[],'ylim',[ydim1 ydim2],'xlim',[xdim1 xdim2],'xtick',xt,'ytick',yt)
box('off')
set(gca,'fontsize',fsize)

sp3 = subplot(713);
set(sp3,'position',[x0 y0+(nfig-3)*(gapy+dy) dx dy])
[i j]=size(d5.FabbriOutput);
p3=plot(d5.FabbriOutput(j).time,d5.FabbriOutput(j).volt);
set(p3,'linewidth',lwid)
set(gca,'xticklabel',[],'ylim',[ydim1 ydim2],'xlim',[xdim1 xdim2],'xtick',xt,'ytick',yt)
box('off')
set(gca,'fontsize',fsize)

sp4 = subplot(714);
set(sp4,'position',[x0 y0+(nfig-4)*(gapy+dy) dx dy])
[i j]=size(d3.FabbriOutput);
p4=plot(d3.FabbriOutput(j).time,d3.FabbriOutput(j).volt);
set(p4,'linewidth',lwid)
set(gca,'xticklabel',[],'ylim',[ydim1 ydim2],'xlim',[xdim1 xdim2],'xtick',xt,'ytick',yt)
box('off')
set(gca,'fontsize',fsize)
ylabel('Voltage (mV)')

sp5 = subplot(715);
set(sp5,'position',[x0 y0+(nfig-5)*(gapy+dy) dx dy])
[i j]=size(d6.FabbriOutput);
p5=plot(d6.FabbriOutput(j).time,d6.FabbriOutput(j).volt);
set(p5,'linewidth',lwid)
set(gca,'xticklabel',[],'ylim',[ydim1 ydim2],'xlim',[xdim1 xdim2],'xtick',xt,'ytick',yt)
box('off')
set(gca,'fontsize',fsize)

sp6 = subplot(716);
set(sp6,'position',[x0 y0+(nfig-6)*(gapy+dy) dx dy])
[i j]=size(d4.FabbriOutput);
p6=plot(d4.FabbriOutput(j).time,d4.FabbriOutput(j).volt);
set(p6,'linewidth',lwid)
set(gca,'xticklabel',[],'ylim',[ydim1 ydim2],'xlim',[xdim1 xdim2],'xtick',xt,'ytick',yt)
box('off')
set(gca,'fontsize',fsize)

sp7 = subplot(717);
set(sp7,'position',[x0 y0+(nfig-7)*(gapy+dy) dx dy])
[i j]=size(d7.FabbriOutput);
p7=plot(d7.FabbriOutput(j).time,d7.FabbriOutput(j).volt);
set(p7,'linewidth',lwid)
set(gca,'ylim',[ydim1 ydim2],'xlim',[xdim1 xdim2],'xtick',xt,'ytick',yt)
box('off')
set(gca,'fontsize',fsize)

xlabel('Time (s)')

set(gca,'fontsize',14)

% print -dtiff -r300dpi APTime2.tiff
