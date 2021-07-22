function biomarkers = Severi_findfiducial(time,volt)
% close all

[peaks, peakpos]=findpeaks(volt);

% figure
% plot(time,volt)

% hold on
% 
% plot(time(peakpos),peaks,'rx')

CL = time(peakpos(end-1))-time(peakpos(end-2));

biomarkers.CL = CL;