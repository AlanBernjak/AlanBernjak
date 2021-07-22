% Script reproduces Figure3b
% Corresponding data are in files 'Figure3b_data'

% All input parameters values of the model are saved in the 'parameters'
% Matlab structure. 

% All features of the AP are saved in the 'biomarkers' structure. 

% All AP time series are saved in the 'FabbriOutput' structure.


clear
close all

load('Figure3b_data')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure; 
hold on

plotno = 1:length(biomarkers);

allgkr = [];

for count = 1:2:length(plotno)
    
    allgkr = [allgkr; parameters(count).x_g_Kr];
        
    ind = plotno(count);
    
    time = FabbriOutput(ind).time;
    volt = FabbriOutput(ind).volt;
    
    %%%
 
    [peaksy, peaksx_ind] = findpeaks(volt,'minpeakheight',0);
    
    if length(peaksx_ind)<2
        return
    end
    
    selpeak = peaksx_ind(end-4);
    
    [peakmin_1, peakminind_1] = min(volt(peaksx_ind(end-5):peaksx_ind(end-4)));
    peakminind_1 = peakminind_1 + peaksx_ind(end-5) -1;
    
    [maxv, maxi] = max(diff(volt(peakminind_1:selpeak)));
    
    maxi = maxi+peakminind_1-1;
    
    %%%
    
    CL_sampl = peaksx_ind(end-1)-peaksx_ind(end-2);

    int2 = CL_sampl*4; 

    pl=plot(time(maxi:maxi+round(int2))-time(maxi),volt(maxi:maxi+round(int2)));
    
    set(pl,'linewidth',1.5)
        
end

legend({'{\it g}_{Kr} = 1.0','{\it g}_{Kr} = 0.9','{\it g}_{Kr} = 0.8','{\it g}_{Kr} = 0.7'})

set(gca,'xlim',[0 1.32],'ylim',[-67 32])
set(gca,'ytick',[-60 -40 -20 0 20],'xtick',[0 0.2 0.4 0.6 0.8 1.0 1.2])

xlabel('Time (s)')
ylabel('Voltage (mV)');

set(gca,'fontsize',14)

set(gca,'linewidth',1)

chH = get(gca,'Children');
set(chH(end),'linewidth',3)

% set(f,'position',[353 76 997 254]);
% t = uitable(f,'Data',rhomat,'Units','Normalized','RowName',fields,'ColumnName',fields,'Position',[0 0 1 1]);
        
% print -dtiff APgKr2x.tiff
        
    
    







    