% Script reproduces Figure8a
% Corresponding data are in files 'Figure8a_data1' to 'Figure8a_data5',

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% All AP time series are saved in the 'FabbriOutput' structure.

function Figure8a_mfile

close all
f1=figure;

allcolor = get(gca,'colororder');

BL=load('Figure8a_data1');

BLsymp=load('Figure8a_data2');
BLsymp_lowg=load('Figure8a_data3');
BLsymp_hypoKo35=load('Figure8a_data4');
BLsymp_hypoKo3=load('Figure8a_data5');

[x, y, maxiBL, CLBL] = findmaxdvdt(BL.FabbriOutput.volt);

[x, y, maxiBLsymp, CLBLsymp] = findmaxdvdt(BLsymp.FabbriOutput.volt);
[x, y, maxiBLsymp_lowg, CLBLsymp_lowg] = findmaxdvdt(BLsymp_lowg.FabbriOutput.volt);
[x, y, maxiBLsymp_hypoKo35, CLBLsymp_hypoKo35] = findmaxdvdt(BLsymp_hypoKo35.FabbriOutput.volt);
[x, y, maxiBLsymp_hypoKo3, CLBLsymp_hypoKo3] = findmaxdvdt(BLsymp_hypoKo3.FabbriOutput.volt);

%
fact = 1.8;

pbl = plot(BL.FabbriOutput.time(maxiBL:maxiBL+round(CLBL*fact))-BL.FabbriOutput.time(maxiBL),BL.FabbriOutput.volt(maxiBL:maxiBL+round(CLBL*fact)));
hold on

ps = plot(BLsymp.FabbriOutput.time(maxiBLsymp:maxiBLsymp+round(CLBLsymp*fact))-BLsymp.FabbriOutput.time(maxiBLsymp),BLsymp.FabbriOutput.volt(maxiBLsymp:maxiBLsymp+round(CLBLsymp*fact)));
pslg = plot(BLsymp_lowg.FabbriOutput.time(maxiBLsymp_lowg:maxiBLsymp_lowg+round(CLBLsymp_lowg*fact))-BLsymp_lowg.FabbriOutput.time(maxiBLsymp_lowg),BLsymp_lowg.FabbriOutput.volt(maxiBLsymp_lowg:maxiBLsymp_lowg+round(CLBLsymp_lowg*fact)));
psh35 = plot(BLsymp_hypoKo35.FabbriOutput.time(maxiBLsymp_hypoKo35:maxiBLsymp_hypoKo35+round(CLBLsymp_hypoKo35*fact))-BLsymp_hypoKo35.FabbriOutput.time(maxiBLsymp_hypoKo35),BLsymp_hypoKo35.FabbriOutput.volt(maxiBLsymp_hypoKo35:maxiBLsymp_hypoKo35+round(CLBLsymp_hypoKo35*fact)));
psh3 = plot(BLsymp_hypoKo3.FabbriOutput.time(maxiBLsymp_hypoKo3:maxiBLsymp_hypoKo3+round(CLBLsymp_hypoKo3*fact))-BLsymp_hypoKo3.FabbriOutput.time(maxiBLsymp_hypoKo3),BLsymp_hypoKo3.FabbriOutput.volt(maxiBLsymp_hypoKo3:maxiBLsymp_hypoKo3+round(CLBLsymp_hypoKo3*fact)));

set(pbl,'linewidth',2,'color',[0.4 0.4 0.4],'linestyle','--')
set(ps,'linewidth',3,'color',allcolor(1,:))
set(pslg,'linewidth',2,'color',allcolor(2,:))
set(psh35,'linewidth',2,'color',allcolor(3,:))
set(psh3,'linewidth',2,'color',allcolor(4,:))


ll=legend('NG, {\itK}o=4',...
          'NG, {\itK}o=4 + Symp.Act.',...
          'LG, {\itK}o=4 + Symp.Act.',...
          'LG, {\itK}o=3.5 + Symp.Act.',...
          'LG, {\itK}o=3 + Symp.Act.');...


set(gcf,'position',[450 350 949 420])

set(gca,'position',[0.09 0.14 0.88 0.78])

set(gca,'box','off')

set(gca,'xlim',[0 2.05],'ylim',[-72 32])

ylabel('Voltage (mV)')
xlabel('Time (s)')

set(gca,'fontsize',18)
set(ll,'fontsize',16,'position',[0.65 0.57 0.3 0.34])

set(gca,'position',[0.09 0.19 0.88 0.73])
set(gcf,'position',[450 350 1094 420])

% print -dtiff -r300 APSymp_x.tiff

end


%%%
 
function [peaksy, peaksx_ind, maxi, CL] = findmaxdvdt(volt)

[peaksy, peaksx_ind] = findpeaks(volt,'minpeakheight',0);

CL=mean(diff(peaksx_ind));

if length(peaksx_ind)<2
    return
end

selpeak = peaksx_ind(end-5);

[peakmin_1, peakminind_1] = min(volt(peaksx_ind(end-6):peaksx_ind(end-5)));
peakminind_1 = peakminind_1 + peaksx_ind(end-6) -1;

[maxv, maxi] = max(diff(volt(peakminind_1:selpeak)));

maxi = maxi+peakminind_1-1;
    
end
