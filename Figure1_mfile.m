% The code reproduces Figure1: Stability of AP feature
% The corresponding data are in the 'Figure1_data.mat' file

clear
close all

load('Figure1_data.mat')

% timecut = 
time = FabbriOutput.time;
volt = FabbriOutput.volt;

time = time(1:75000);
volt = volt(1:75000);

figure
plot(time,volt,'bx')

time2 = linspace(time(1),time(end),10000000);
volt2 = spline(time,volt,time2);

time = time2;
volt = volt2;

hold on
plot(time2,volt2,'r')

% 
% figure
% 
% plot(time2,volt2,'bx');
% hold on
% plot(time,volt,'r')


[maxpks, maxlocs] = findpeaks(volt,'MinPeakHeight',1);

hold on
plot(time(maxlocs),maxpks,'rx')

% figure
% plot(time,volt*(-1))

[minpks, minlocs] = findpeaks(volt*(-1),'MinPeakHeight',50);
minpks = volt(minlocs);

plot(time(minlocs),minpks,'gx')

dvolt = diff(volt);

figure
plot(diff(volt))

[maxdiff, maxdifflocs] = findpeaks(dvolt,'MinPeakDistance',100000);

hold on

plot(maxdifflocs,maxdiff,'rx')

%%%%%%%%%%%%%%%

dmaxpeaks = diff(maxpks);
dminpeaks = diff(minpks);

dtmaxlocs = diff(time(maxlocs));

dtminlocs = diff(time(minlocs));

dCL = diff(dtmaxlocs);

% maxpks
% 
% minpks

%%%%

f1 = figure;

p1=plot(time(maxlocs(1:end-2)), dCL./median(dtmaxlocs));
hold on
p2=plot(time(maxlocs(1:end-1)), dmaxpeaks./median(maxpks));
p3=plot(time(minlocs(1:end-1)), dminpeaks./median(minpks));

% ylim = [-0.0002 0.0005];

set(gca,'xlim',[-2 62],'ylim',ylim)

set(gca,'linewidth',1)

set([p1 p2 p3],'linewidth',2)

yl = ylabel('Relative beat-to-beat difference');
xl = xlabel('Time (s)');

xlp = get(xl,'position');
% set(xl,'position',[120 -0.0018 xlp(3)])

set(gca,'fontsize',14)
set(gca,'box','off')

% sp2=subplot(122);
% 
% p3=plot(time(maxlocs(1:end-2)), dCL./median(dtmaxlocs));
% hold on
% p1=plot(time(maxlocs(1:end-1)), dmaxpeaks./median(maxpks));
% p2=plot(time(minlocs(1:end-1)), dminpeaks./median(minpks));
% 
% set(gca,'xlim',[500 600],'ylim',ylim,'yticklabel',[])
% set(gca,'linewidth',1)
% 
% set([p1 p2 p3],'linewidth',2)

legend({'\DeltaCL/median(CL)','\DeltaPP/median(PP)','\DeltaMDP/median(MDP)'})

set(gca,'fontsize',14)

dx=0.02;
ypos = 0.15;
yh = 0.77;

pos1 = get(gca,'position');

set(gca,'position',[pos1(1) 0.185 0.8 0.72]);
% set(sp2,'position',[pos2(1)-dx ypos pos2(3)+dx yh]);

% set(gcf,'position',[680 558 525 345])

set(gca,'box','off')
% [mean(dmaxpeaks(20:60)) std(dmaxpeaks(20:60));
%  mean(dmaxpeaks(310:350)) std(dmaxpeaks(310:350));
%  mean(dmaxpeaks(610:650)) std(dmaxpeaks(610:650))];
    
set(gca,'fontsize',14)

set(gca,'linewidth',1)

% print -dtiff Stability3x.tiff




