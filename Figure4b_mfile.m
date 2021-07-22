% Script reproduces Figure4b
% Corresponding data are in files 'Figure4b_data'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% All AP time series are saved in the 'FabbriOutput' structure.

clear
close all

load('Figure4b_data')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure; 
hold on

plotno = 1:2:length(biomarkers);

allv12 = 0;

for count = 1:length(plotno)
    
    allv12 = [allv12; parameters(count).IHerg_shift_pa];
    
    ind = plotno(count);
    
    time = FabbriOutput(ind).time;
    volt = FabbriOutput(ind).volt;
    
    %%%
 
    [peaksy, peaksx_ind] = findpeaks(volt,'minpeakheight',0);
    
    if length(peaksx_ind)<2
        return
    end
    
    selpeak = peaksx_ind(end-5);
    
    [peakmin_1, peakminind_1] = min(volt(peaksx_ind(end-6):peaksx_ind(end-5)));
    peakminind_1 = peakminind_1 + peaksx_ind(end-6) -1;
    
    [maxv, maxi] = max(diff(volt(peakminind_1:selpeak)));
    
    maxi = maxi+peakminind_1-1;
    
    %%%
    
    CL_sampl = peaksx_ind(end-1)-peaksx_ind(end-2);

    int2 = CL_sampl*5; 

    pl=plot(time(maxi:maxi+round(int2))-time(maxi),volt(maxi:maxi+round(int2)));
    
    set(pl,'linewidth',1.5)
        
end

legpos = legend({'IVshift = 0 mV','IVshift = -1.0 mV','IVshift = -2.0 mV',...
    'IVshift = -3.0 mV'});

set(legpos,'position',[0.1905    0.6933    0.3179    0.2226])

set(gca,'xlim',[0 1.82],'ylim',[-81 32])
set(gca,'ytick',[-80 -60 -40 -20 0 20],'xtick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8])

xlabel('Time (s)')
ylabel('Voltage (mV)');

set(gca,'fontsize',14)

set(gca,'linewidth',1)

chH = get(gca,'Children');
set(chH(end),'linewidth',3)

% set(gca,'Children',[chH(end);chH(1:end-1)])
    
% print -dtiff APV12x2.tiff
        
    
    







    