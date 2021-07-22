% plots the AP waveform and AP features saved in the 'testFabbri' structure

clear
close all

load testFabbri

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure; set(f1,'position',[84 188 560 420])
hold on
% f2 = figure; set(f2,'position',[646 189 560 420])
% hold on

plotno = 1:length(biomarkers);

for count = 1:length(plotno)
    
    ind = plotno(count);
    
    time = FabbriOutput(ind).time;
    volt = FabbriOutput(ind).volt;
    
    figure(f1)
    plot(time,volt)
    
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
    
%     int1 = CL_sampl*0.5;
    int2 = CL_sampl*5; 
%     int2 = 2000;
    
%     figure(f2)
%     
%     plot(time(maxi:maxi+round(int2))-time(maxi),volt(maxi:maxi+round(int2)))
    
    %%%%
%     for stleg = 1 : length(fieldnames(parameters))-1
%         fnames = fieldnames(parameters);
%         leglabel{stleg} = fnames{stleg};
%         fval = parameters{stleg};
%         
%     end
%     
%     
end

biom_temp = findfiducial4(time,volt,0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

