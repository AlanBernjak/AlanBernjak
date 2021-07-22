% Script reproduces Figure6bcd
% Corresponding data are in files 'Figure6bcd_data1','Figure6bcd_data2',
% 'Figure6bcd_data3' and 'Figure6bcd_data4'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% All AP time series are saved in the 'FabbriOutput' structure.


clear
close all

f1 = figure;

xlength = 0.84;
yheight = 0.27;

set(f1,'position',[520 100 560 630])

sp1=subplot(311);
hold on
set(sp1,'position',[0.13 0.71 xlength yheight]);

d1=load('Figure6bcd_data1');
load('Figure6bcd_data2');

Output(1)=d1.FabbriOutput;
Output(2)=FabbriOutput(1);

plotno = length(Output);

for count = 1:plotno
    
%     ind = allK(count);
    
    time = Output(count).time;
    volt = Output(count).volt;
    
    %%%
 
    [peaksy, peaksx_ind] = findpeaks(volt,'minpeakheight',0);
    
    if length(peaksx_ind)<2
        return
    end
    
    selpeak = peaksx_ind(end-3);
    
    [peakmin_1, peakminind_1] = min(volt(peaksx_ind(end-4):peaksx_ind(end-3)));
    peakminind_1 = peakminind_1 + peaksx_ind(end-4) -1;
    
    [maxv, maxi] = max(diff(volt(peakminind_1:selpeak)));
    
    maxi = maxi+peakminind_1-1;
    
    %%%
    
    CL_sampl = peaksx_ind(end-1)-peaksx_ind(end-2);

    int2 = CL_sampl*2; 
    
    pl=plot(time(maxi:maxi+round(int2))-time(maxi),volt(maxi:maxi+round(int2)));
    
    set(pl,'linewidth',1.5)
    
    clear maxv maxi peakmin* peaksy peaks_ind time volt
        
end

ll1=legend({'Normoglycaemia','Low glucose'});

set(ll1,'position',[0.28 0.89 0.3 0.09],'fontsize',12)

set(gca,'xlim',[-0.01 1.42],'ylim',[-72 32])
set(gca,'ytick',[-60 -40 -20 0 20],'xtick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4],'xticklabel',[])

% xlabel('Time (s)')
ylabel('Voltage (mV)');

set(gca,'fontsize',14)

set(gca,'linewidth',1)
    
% print -dtiff -r300dpi KoAP3.tiff
        
clear biomarkers FabbriOutput parameters

%%

sp2=subplot(312);
set(sp2,'position',[0.13 0.41 xlength yheight]);
hold on

d1=load('Figure6bcd_data1');
load('Figure6bcd_data3');

Output(1)=FabbriOutput(1);
Output(2)=FabbriOutput(2);
Output(3)=FabbriOutput(3);

plotno = length(Output);

for count = 1:plotno
    
%     ind = allK(count);
    
    time = Output(count).time;
    volt = Output(count).volt;
    
    %%%
 
    [peaksy, peaksx_ind] = findpeaks(volt,'minpeakheight',0);
    
    if length(peaksx_ind)<2
        return
    end
    
    selpeak = peaksx_ind(end-3);
    
    [peakmin_1, peakminind_1] = min(volt(peaksx_ind(end-4):peaksx_ind(end-3)));
    peakminind_1 = peakminind_1 + peaksx_ind(end-4) -1;
    
    [maxv, maxi] = max(diff(volt(peakminind_1:selpeak)));
    
    maxi = maxi+peakminind_1-1;
    
    %%%
    
    CL_sampl = peaksx_ind(end-1)-peaksx_ind(end-2);

    int2 = CL_sampl*2; 
    
    pl=plot(time(maxi:maxi+round(int2))-time(maxi),volt(maxi:maxi+round(int2)));
    
    set(pl,'linewidth',1.5)
    
    clear maxv maxi peakmin* peaksy peaks_ind time volt
        
end

ll2=legend({'{\it K}o = 4.0 mmol/L','{\it K}o = 3.5 mmol/L','{\it K}o = 3.0 mmol/L'});

set(ll2,'position',[0.28 0.538 0.3 0.09],'fontsize',12)

t12=text(0.23,29,'Normoglycaemia');

set(gca,'xlim',[-0.01 1.42],'ylim',[-72 32])
set(gca,'ytick',[-60 -40 -20 0 20],'xtick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4],'xticklabel',[])

% xlabel('Time (s)')
ylabel('Voltage (mV)');

set(gca,'fontsize',14)
set(t12,'fontsize',14)

set(gca,'linewidth',1)
        
clear biomarkers FabbriOutput parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%

sp3 = subplot(313);
set(sp3,'position',[0.13 0.11 xlength yheight]);
hold on

load('Figure6bcd_data4');

Output(1)=FabbriOutput(1);
Output(2)=FabbriOutput(2);
Output(3)=FabbriOutput(3);

plotno = length(Output);

for count = 1:plotno
    
%     ind = allK(count);
    
    time = Output(count).time;
    volt = Output(count).volt;
    
    %%%
 
    [peaksy, peaksx_ind] = findpeaks(volt,'minpeakheight',0);
    
    if length(peaksx_ind)<2
        return
    end
    
    selpeak = peaksx_ind(end-3);
    
    [peakmin_1, peakminind_1] = min(volt(peaksx_ind(end-4):peaksx_ind(end-3)));
    peakminind_1 = peakminind_1 + peaksx_ind(end-4) -1;
    
    [maxv, maxi] = max(diff(volt(peakminind_1:selpeak)));
    
    maxi = maxi+peakminind_1-1;
    
    %%%
    
    CL_sampl = peaksx_ind(end-1)-peaksx_ind(end-2);

    int2 = CL_sampl*2; 
    
    pl=plot(time(maxi:maxi+round(int2))-time(maxi),volt(maxi:maxi+round(int2)));
    
    set(pl,'linewidth',1.5)
    
    clear maxv maxi peakmin* peaksy peaks_ind time volt
        
end

ll3=legend({'{\it K}o = 4.0 mmol/L','{\it K}o = 3.5 mmol/L','{\it K}o = 3.0 mmol/L'});

set(ll3,'position',[0.28 0.237 0.3 0.09],'fontsize',12)

t13=text(0.23,29,'Low glucose');

set(gca,'xlim',[-0.01 1.42],'ylim',[-72 32])
set(gca,'ytick',[-60 -40 -20 0 20],'xtick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4],'xticklabel',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4])

xlabel('Time (s)')
ylabel('Voltage (mV)');

set(gca,'fontsize',14)
set(t13,'fontsize',14)

set(gca,'linewidth',1)

% print -dtiff -r300dpi AP_KoLowGKo2.tiff






% % plot(d1.FabbriOutput.time,d1.FabbriOutput.volt)
% % plot(d2.FabbriOutput.time,d2.FabbriOutput.volt)
% 
% 
% plotno = length(FabbriOutput);
% 
% for count = 1:plotno
%     
% %     ind = allK(count);
%     
%     time = FabbriOutput(count).FabbriOutput.time;
%     volt = FabbriOutput(count).FabbriOutput.volt;
%     
%     %%%
%  
%     [peaksy, peaksx_ind] = findpeaks(volt,'minpeakheight',0);
%     
%     if length(peaksx_ind)<2
%         return
%     end
%     
%     selpeak = peaksx_ind(end-3);
%     
%     [peakmin_1, peakminind_1] = min(volt(peaksx_ind(end-4):peaksx_ind(end-3)));
%     peakminind_1 = peakminind_1 + peaksx_ind(end-4) -1;
%     
%     [maxv, maxi] = max(diff(volt(peakminind_1:selpeak)));
%     
%     maxi = maxi+peakminind_1-1;
%     
%     %%%
%     
%     CL_sampl = peaksx_ind(end-1)-peaksx_ind(end-2);
% 
%     int2 = CL_sampl*2; 
%     
%     pl=plot(time(maxi:maxi+round(int2))-time(maxi),volt(maxi:maxi+round(int2)));
%     
%     set(pl,'linewidth',1.5)
%         
% end
% 
% % ll=legend({'Normoglycaemia, Ko = 4.0',...
% %     'Low glucose & Ko = 4.0',...
% %     'Low glucose & Ko = 3.5',...
% %     'Low glucose & Ko = 3.0'});
% 
% chH = get(gca,'Children');
% set(chH(4),'linewidth',3)
% 
% set(gca,'xlim',[0 1.3],'ylim',[-68 33])
% set(gca,'ytick',[-60 -40 -20 0 20],'xtick',[0 0.2 0.4 0.6 0.8 1.0 1.2])
% 
% xlabel('Time (s)')
% ylabel('Voltage (mV)')
% 
% ll2=legend({'NG, {\it K}o = 4.0 mmol/L','LG, {\it K}o = 4.0 mmol/L','LG, {\it K}o = 3.5 mmol/L','LG, {\it K}o = 3.0 mmol/L'});
% 
% set(ll2,'position',[0.2655 0.37 0.3232 0.1690])
% 
% % t2=text(0.2,25,'Low glucose');
% % 
% set([t1],'fontsize',14)
% 
% % set(ll,'fontsize',12,'position',[0.21 0.71 0.42 0.21])
% set(gca,'fontsize',14)
% 
% 
% % print -dtiff -r300dpi AP_KoLowGKo.tiff

    
    







    